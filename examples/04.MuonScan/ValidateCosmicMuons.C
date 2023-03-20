#include <TRestGeant4Event.h>

Int_t ValidateCosmicMuons(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    TRestGeant4Metadata* md = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");

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
    constexpr double averageTotalEnergyRef = 15.244;

    double averageSensitiveEnergy = 0;
    constexpr double averageSensitiveEnergyRef = 14.8731;

    constexpr double tolerance = 0.1;

    Int_t detDown = 0;
    Int_t detUp = 0;
    Int_t detBoth = 0;
    for (Long64_t i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);

        averageTotalEnergy += event->GetTotalDepositedEnergy() / nEvents;
        averageSensitiveEnergy += event->GetSensitiveVolumeEnergy() / nEvents;

        Bool_t down = !TMath::IsNaN(event->GetFirstPositionInVolume(md->GetActiveVolumeID("det_dw_01")).X());
        Bool_t up = !TMath::IsNaN(event->GetFirstPositionInVolume(md->GetActiveVolumeID("det_up_01")).X());
        if (down) detDown++;
        if (up) detUp++;
        if (down && up) detBoth++;
    }

    cout << "Number of events that crossed down detector : " << detDown << endl;
    cout << "Number of events that crossed up detector : " << detUp << endl;
    cout << "Number of events that crossed both detectors : " << detBoth << endl;

    if (detDown != 10000) {
        cout << "The number of cosmics crossing the sensitive volume is not 10000!" << endl;
        return 13;
    }

    if (detUp != 462) {
        cout << "The number of cosmics crossing the up volume is not 462!" << endl;
        return 14;
    }

    if (detBoth != 462) {
        cout << "The number of cosmics crossing both volumes is not 462!" << endl;
        return 15;
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
