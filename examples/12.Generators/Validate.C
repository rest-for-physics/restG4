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
    if (run.GetEntries() != 200000) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 2;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    TVector3 sourceDirection = metadata->GetParticleSource()->GetDirection();
    cout << "original source direction: (" << sourceDirection.X() << ", " << sourceDirection.Y() << ", "
         << sourceDirection.Z() << ")" << endl;

    TH1D thetaHist("thetaHist", "Theta angle from source direction", 100, 0, TMath::Pi());
    for (int i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);
        TVector3 direction = event->GetPrimaryEventDirection();
        const auto theta = sourceDirection.Angle(direction);
        thetaHist.Fill(theta);
    }
    thetaHist.Scale(1. / thetaHist.Integral(), "width");

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
        return 3;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
