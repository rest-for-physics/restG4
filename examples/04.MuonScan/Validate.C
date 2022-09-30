
#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestRun.h>

Int_t Validate(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);

    if (run.GetParentRunNumber() != 0) {
        cout << "Parent run number value : " << run.GetParentRunNumber() << endl;
        cout << "The parent run number from restG4 generated file should be 0" << endl;
        return 1;
    }

    if (run.GetRunNumber() != 111) {
        cout << "Run number value : " << run.GetRunNumber() << endl;
        cout << "The run number on the Muon example should be 111 by default!" << endl;
        return 2;
    }

    if (run.GetRunType() != "restG4") {
        cout << "Run type : " << run.GetRunType() << endl;
        cout << "The run type of restG4 generated data should be 'restG4'!" << endl;
        return 3;
    }

    if (run.GetRunTag() != "MuonsFromPoint") {
        cout << "Run tag : " << run.GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'MuonsFromPoint'!" << endl;
        return 4;
    }

    if (run.GetEntries() != 1000) {
        cout << "Run entries : " << run.GetEntries() << endl;
        cout << "The number of stored events should match the reference value of 1000" << endl;
        return 5;
    }

    cout << "Testing reading of Geant4 metadata class" << endl;
    TRestGeant4Metadata* geant4Metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    if (!geant4Metadata) {
        cout << "Problem reading Geant4 metadata class!" << endl;
        return 6;
    }
    geant4Metadata->PrintMetadata();

    if (geant4Metadata->GetNumberOfActiveVolumes() != 2) {
        cout << "The number of active volumes does not match the reference value of 2" << endl;
        return 7;
    }

    if (geant4Metadata->GetMaterialsReference() != "0.0") {
        cout << "Materials version was not defined!" << endl;
        cout << "The parent run number from restG4 generated file should be 0" << endl;
        return 8;
    }

    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    double nEvents = run.GetEntries();

    double averageTotalEnergy = 0;
    constexpr double averageTotalEnergyRef = 18.074;

    double averageSensitiveEnergy = 0;
    constexpr double averageSensitiveEnergyRef = 8.96635;

    double averageNumberOfTracks = 0;
    constexpr double averageNumberOfTracksRef = 393.195;

    double averageNumberOfHitsVolume0 = 0;
    constexpr double averageNumberOfHitsVolume0Ref = 66.076;

    double averageNumberOfHitsVolume1 = 0;
    constexpr double averageNumberOfHitsVolume1Ref = 59.981;

    TVector3 averagePosition = {};
    const TVector3 averagePositionRef = {-0.59902, 0.0256995, 300.385};

    constexpr double tolerance = 0.001;

    for (Long64_t i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);

        averageTotalEnergy += event->GetTotalDepositedEnergy() / nEvents;
        averageSensitiveEnergy += event->GetSensitiveVolumeEnergy() / nEvents;
        averageNumberOfHitsVolume0 += event->GetNumberOfHits(0) / nEvents;
        averageNumberOfHitsVolume1 += event->GetNumberOfHits(1) / nEvents;
        averageNumberOfTracks += event->GetNumberOfTracks() / nEvents;
        averagePosition += event->GetMeanPositionInVolume(0) * (1.0 / nEvents);
    }

    cout << "Average total energy: " << averageTotalEnergy << " keV" << endl;
    cout << "Average sensitive energy: " << averageSensitiveEnergy << " keV" << endl;
    cout << "Average number of hits in volume 0: " << averageNumberOfHitsVolume0 << endl;
    cout << "Average number of hits in volume 1: " << averageNumberOfHitsVolume1 << endl;
    cout << "Average number of tracks: " << averageNumberOfTracks << endl;
    cout << "Average position: (" << averagePosition.x() << ", " << averagePosition.y() << ", "
         << averagePosition.z() << ") mm" << endl;

    if (TMath::Abs(averageNumberOfTracks - averageNumberOfTracksRef) / averageNumberOfTracksRef > tolerance) {
        cout << "The average number of tracks does not match the reference value of "
             << averageNumberOfTracksRef << endl;
        return 9;
    }

    if (TMath::Abs(averageNumberOfHitsVolume0 - averageNumberOfHitsVolume0Ref) /
            averageNumberOfHitsVolume0Ref >
        tolerance) {
        cout << "The average number of hits in volume 0 does not match the reference value of "
             << averageNumberOfHitsVolume0Ref << endl;
        return 10;
    }

    if (TMath::Abs(averageNumberOfHitsVolume1 - averageNumberOfHitsVolume1Ref) /
            averageNumberOfHitsVolume1Ref >
        tolerance) {
        cout << "The average number of hits in volume 1 does not match the reference value of "
             << averageNumberOfHitsVolume1Ref << endl;
        return 11;
    }

    if (TMath::Abs(averageSensitiveEnergy - averageSensitiveEnergyRef) / averageSensitiveEnergyRef >
        tolerance) {
        cout << "The average sensitive volume energy does not match the reference value of "
             << averageSensitiveEnergyRef << endl;
        return 12;
    }

    if (TMath::Abs(averageTotalEnergy - averageTotalEnergyRef) / averageTotalEnergyRef > tolerance) {
        cout << "The average total energy does not match the reference value of " << averageTotalEnergyRef
             << endl;
        return 13;
    }

    if (TMath::Abs(averagePosition.Mag() - averagePositionRef.Mag()) / averagePositionRef.Mag() > tolerance) {
        cout << "The average position does not match the reference value of "
             << "(" << averagePositionRef.x() << ", " << averagePositionRef.y() << ", "
             << averagePositionRef.z() << ") mm" << endl;
        return 14;
    }

    return 0;
}
