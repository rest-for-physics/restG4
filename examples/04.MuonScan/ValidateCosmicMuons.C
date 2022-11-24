#include <TRestGeant4Event.h>

Int_t ValidateCosmicMuons(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    if (run.GetRunTag() != "CosmicMuons") {
        cout << "Run tag: " << run.GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'CosmicMuons'" << endl;
        return 1;
    }

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() != 10000) {
        cout << "The number of entries is wrong!" << endl;
        cout << "Number of entries : " << run.GetEntries() << endl;
        return 2;
    }

    double nEvents = run.GetEntries();

    double averageTotalEnergy = 0;
    constexpr double averageTotalEnergyRef = 14.0729;

    double averageSensitiveEnergy = 0;
    constexpr double averageSensitiveEnergyRef = 13.5387;

    constexpr double tolerance = 0.1;

    for (Long64_t i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);

        averageTotalEnergy += event->GetTotalDepositedEnergy() / nEvents;
        averageSensitiveEnergy += event->GetSensitiveVolumeEnergy() / nEvents;
    }

    cout << "Average total energy: " << averageTotalEnergy << " keV" << endl;
    cout << "Average sensitive energy: " << averageSensitiveEnergy << " keV" << endl;

    if (TMath::Abs(averageSensitiveEnergy - averageSensitiveEnergyRef) / averageSensitiveEnergyRef >
        tolerance) {
        cout << "The average sensitive volume energy does not match the reference value of "
             << averageSensitiveEnergyRef << endl;
        return 3;
    }

    if (TMath::Abs(averageTotalEnergy - averageTotalEnergyRef) / averageTotalEnergyRef > tolerance) {
        cout << "The average total energy does not match the reference value of " << averageTotalEnergyRef
             << endl;
        return 4;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";

    return 0;
}
