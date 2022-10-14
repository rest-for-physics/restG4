#include <TRestGeant4Event.h>

Int_t Validate(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);

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

    TH1D thetaHist("thetaHist", "Theta angle from source direction", 100, 0, TMath::Pi());
    TH1D energyHist("energyHist", "Primary muon energy", 200, 0, 1E9);
    for (int i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);
        TVector3 direction = event->GetPrimaryEventDirection();
        const auto theta = sourceDirection.Angle(direction);
        thetaHist.Fill(theta);

        Double_t energy = event->GetPrimaryEventEnergy();
        energyHist.Fill(energy);
    }
    thetaHist.Scale(1. / thetaHist.Integral(), "width");
    energyHist.Scale(1. / energyHist.Integral(), "width");

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
