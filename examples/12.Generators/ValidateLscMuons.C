#include <TRestGeant4Event.h>

Int_t ValidateLscMuons(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);

    cout << "Run tag: " << run.GetRunTag() << endl;
    if (run.GetRunTag() != "LscMuons") {
        cout << "The run tag of the basic validation test should be 'LscMuons'!" << endl;
        return 1;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    const Int_t volumeId = metadata->GetActiveVolumeID("gas");
    if (volumeId < 0) {
        cout << "The active volume 'gas' was not found!" << endl;
        return 2;
    }

    cout << "Generated primaries: " << metadata->GetNumberOfEvents() << endl;
    if (metadata->GetNumberOfEvents() != 20000) {
        cout << "Bad number of generated primaries: " << metadata->GetNumberOfEvents() << endl;
        return 3;
    }

    // The muon crosses the active volume in a straight line, so the track length is the distance between
    // the first and the last hit it leaves inside the volume
    TVector3 sourceDirection = metadata->GetParticleSource()->GetDirection();
    TH1D thetaHist("thetaHist", "Zenith angle from source direction", 50, 0, TMath::Pi());
    vector<double> lengths;
    for (int i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);
        thetaHist.Fill(sourceDirection.Angle(event->GetPrimaryEventDirection()));
        for (unsigned int t = 0; t < event->GetNumberOfTracks(); t++) {
            const TRestGeant4Track& track = event->GetTrack(t);
            if (track.GetParentID() != 0 || track.GetParticleName() != "mu-") {
                continue;
            }
            const TRestGeant4Hits& hits = track.GetHits();
            TVector3 firstHit, lastHit;
            bool found = false;
            for (unsigned int n = 0; n < hits.GetNumberOfHits(); n++) {
                if (hits.GetVolumeId(n) != volumeId) {
                    continue;
                }
                const TVector3 position(hits.GetX(n), hits.GetY(n), hits.GetZ(n));
                if (!found) {
                    firstHit = position;
                    found = true;
                }
                lastHit = position;
            }
            if (found) {
                lengths.push_back((lastHit - firstHit).Mag() / 10.);  // cm
            }
        }
    }
    thetaHist.Scale(1. / thetaHist.Integral(), "width");

    sort(lengths.begin(), lengths.end());
    double lengthAverage = 0;
    for (const auto& length : lengths) {
        lengthAverage += length / lengths.size();
    }
    const double lengthPercentile99 = lengths[(size_t)(0.99 * lengths.size())];

    constexpr double lengthAverageRef = 15.97, lengthPercentile99Ref = 30.28;
    constexpr double tolerance = 0.05;

    cout << "Muons crossing the active volume: " << lengths.size() << endl;
    cout << "Average track length (cm): " << lengthAverage << endl;
    cout << "99% percentile of the track length (cm): " << lengthPercentile99 << endl;
    cout << "Maximum track length (cm): " << lengths.back() << endl;

    if (TMath::Abs(lengthAverage - lengthAverageRef) / lengthAverageRef > tolerance) {
        cout << "The average track length is wrong!" << endl;
        return 4;
    }
    if (TMath::Abs(lengthPercentile99 - lengthPercentile99Ref) / lengthPercentile99Ref > tolerance) {
        cout << "The tail of the track length distribution is wrong!" << endl;
        return 5;
    }

    // The absolute rate is not given by Geant4, it comes from the reference flux at the LSC
    const double days = metadata->GetEquivalentSimulatedTime() / (24. * 60. * 60.);
    const double ratePerDay = lengths.size() / days;
    constexpr double ratePerDayRef = 30.57;

    cout << "Equivalent live time (days): " << days << endl;
    cout << "Rate of muons crossing the active volume (per day): " << ratePerDay << endl;

    if (TMath::Abs(ratePerDay - ratePerDayRef) / ratePerDayRef > tolerance) {
        cout << "The rate of muons crossing the active volume is wrong!" << endl;
        return 6;
    }

    // The sampled zenith distribution must follow the normalized sin(theta) * cos^6(theta)
    auto lscMuonsLambda = [](double* xs, double* ps) {
        if (xs[0] >= 0 && xs[0] <= TMath::Pi() / 2) {
            return 7.0 * TMath::Sin(xs[0]) * TMath::Power(TMath::Cos(xs[0]), 6);
        }
        return 0.0;
    };
    const char* title = "AngularDistribution: LscMuons";
    auto referenceFunction = TF1(title, lscMuonsLambda, 0.0, TMath::Pi());
    referenceFunction.SetTitle(title);

    double diff = 0;
    for (int i = 1; i < thetaHist.GetNbinsX(); i++) {
        diff += TMath::Abs(thetaHist.GetBinContent(i) - referenceFunction.Eval(thetaHist.GetBinCenter(i))) /
                thetaHist.GetNbinsX();
    }
    cout << "Computed difference from reference function: " << diff << endl;
    if (diff > 0.01) {
        cout << "The computed difference from reference function is too large!" << endl;
        return 7;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
