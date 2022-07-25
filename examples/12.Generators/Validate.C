#include <TRestGeant4Event.h>

Int_t Validate(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);

    cout << "Run tag: " << run.GetRunTag() << endl;
    if (run.GetRunTag() != "CosineSquaredCircle") {
        cout << "The run tag of the basic validation test should be 'CosineSquaredCircle'!" << endl;
        return 1;
    }

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() != 100000) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 2;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    TVector3 sourceDirection = metadata->GetParticleSource()->GetDirection();
    cout << "original source direction: (" << sourceDirection.X() << ", " << sourceDirection.Y() << ", "
         << sourceDirection.Z() << ")" << endl;

    double radiusAverage = 0, radiusMin = TMath::Infinity(), radiusMax = 0;
    constexpr double radiusAverageRef = 2.67045, radiusMinRef = 0.0, radiusMaxRef = 4.0;
    constexpr double tolerance = 0.1;

    double energyPrimaryAverage = 0, energyPrimaryMin = TMath::Infinity(), energyPrimaryMax = 0;
    constexpr double energyPrimaryAverageRef = 53818.8, energyPrimaryMinRef = 5000.45,
                     energyPrimaryMaxRef = 149998.0;

    TH1D thetaHist("thetaHist", "Theta angle from source direction", 100, 0, TMath::Pi());
    for (int i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);
        TVector3 direction = event->GetPrimaryEventDirection();
        const auto theta = sourceDirection.Angle(direction);
        thetaHist.Fill(theta);

        Double_t x = event->GetPrimaryEventOrigin().X();
        Double_t z = event->GetPrimaryEventOrigin().Z();
        const auto r = TMath::Sqrt(x * x + z * z) / 100;  // cm
        radiusAverage += r / run.GetEntries();
        if (r < radiusMin) {
            radiusMin = r;
        }
        if (r > radiusMax) {
            radiusMax = r;
        }

        Double_t energy = event->GetPrimaryEventEnergy();
        energyPrimaryAverage += energy / run.GetEntries();
        if (energy < energyPrimaryMin) {
            energyPrimaryMin = energy;
        }
        if (energy > energyPrimaryMax) {
            energyPrimaryMax = energy;
        }
    }
    thetaHist.Scale(1. / thetaHist.Integral(), "width");

    cout << "Average radius (cm): " << radiusAverage << endl;
    cout << "Minimum radius (cm): " << radiusMin << endl;
    cout << "Maximum radius (cm): " << radiusMax << endl;

    cout << "Average energy (keV): " << energyPrimaryAverage << endl;
    cout << "Minimum energy (keV): " << energyPrimaryMin << endl;
    cout << "Maximum energy (keV): " << energyPrimaryMax << endl;

    if (TMath::Abs(radiusAverage - radiusAverageRef) > tolerance) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "radiusAverage (cm): " << radiusAverage << endl;
        return 3;
    }
    if (TMath::Abs(radiusMin - radiusMinRef) > tolerance) {
        cout << "The minimum radius of the distribution is wrong!" << endl;
        cout << "radiusMin (cm): " << radiusMin << endl;
        return 4;
    }
    if (TMath::Abs(radiusMax - radiusMaxRef) > tolerance) {
        cout << "The maximum radius of the distribution is wrong!" << endl;
        cout << "radiusMax (cm): " << radiusMax << endl;
        return 5;
    }

    if (TMath::Abs(energyPrimaryAverage - energyPrimaryAverageRef) > tolerance) {
        cout << "The average energy of the distribution is wrong!" << endl;
        cout << "radiusAverage (cm): " << radiusAverage << endl;
        return 6;
    }
    if (TMath::Abs(energyPrimaryMin - energyPrimaryMinRef) > tolerance) {
        cout << "The minimum energy of the distribution is wrong!" << endl;
        cout << "radiusMin (cm): " << radiusMin << endl;
        return 7;
    }
    if (TMath::Abs(energyPrimaryMax - energyPrimaryMaxRef) > tolerance) {
        cout << "The maximum energy of the distribution is wrong!" << endl;
        cout << "radiusMax (cm): " << radiusMax << endl;
        return 8;
    }

    auto cos2Lambda = [](double* xs, double* ps) {
        if (xs[0] >= 0 && xs[0] <= TMath::Pi() / 2) {
            return TMath::Power(TMath::Cos(xs[0]), 2) / (TMath::Pi() / 4.0);
        }
        return 0.0;
    };
    const char* title = "AngularDistribution: Cos2";
    auto referenceFunction = TF1(title, cos2Lambda, 0.0, TMath::Pi());
    referenceFunction.SetTitle(title);

    TCanvas* canvas = new TCanvas();
    thetaHist.Draw();
    referenceFunction.Draw("same");
    canvas->Print("/tmp/thetaHist.pdf");

    double diff = 0;
    for (int i = 1; i < thetaHist.GetNbinsX(); i++) {
        diff += TMath::Abs(thetaHist.GetBinContent(i) - referenceFunction.Eval(thetaHist.GetBinCenter(i))) /
                thetaHist.GetNbinsX();
    }
    cout << "Computed difference from reference function: " << diff << endl;
    if (diff > 0.01) {
        cout << "The computed difference from reference function is too large!" << endl;
        return 9;
    }

    // Validate energy distribution

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
