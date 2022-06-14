
#include <TRestGeant4Event.h>

Int_t Validate(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);

    if (run.GetParentRunNumber() != 0) {
        cout << "Parent run number value : " << run.GetParentRunNumber() << endl;
        cout << "The parent run number from restG4 generated file should be 0" << endl;
        return 1;
    }

    if (run.GetRunNumber() != 1) {
        cout << "Run number value : " << run.GetRunNumber() << endl;
        cout << "The run number on the validation chain should be 1 by default!" << endl;
        return 2;
    }

    if (run.GetRunType() != "restG4") {
        cout << "Run type : " << run.GetRunType() << endl;
        cout << "The run type of restG4 generated data should be 'restG4'!" << endl;
        return 3;
    }

    if (run.GetRunTag() != "NLDBD") {
        cout << "Run tag : " << run.GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'NLDBD'!" << endl;
        return 4;
    }

    if (run.GetEntries() != 100) {
        cout << "Run entries : " << run.GetEntries() << endl;
        cout << "The number of stored events should match the reference value of 100" << endl;
        return 5;
    }

    cout << "Testing reading of Geant4 metadata class" << endl;
    TRestGeant4Metadata* geant4Metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    if (!geant4Metadata) {
        cout << "Problem reading Geant4 metadata class!" << endl;
        return 6;
    }
    geant4Metadata->PrintMetadata();

    const bool isReferenceGeant4Version = geant4Metadata->GetGeant4Version() == "10.4.3";

    if (geant4Metadata->GetNumberOfActiveVolumes() != 1) {
        cout << "The number of registered does not match the reference value of 1" << endl;
        return 7;
    }

    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    double nEvents = run.GetEntries();

    double averageTotalEnergy = 0;
    const double averageTotalEnergyRef = (!isReferenceGeant4Version) ? 2280.96 : 2221.55;

    double averageSensitiveEnergy = 0;
    const double averageSensitiveEnergyRef = (!isReferenceGeant4Version) ? 2280.96 : 2221.55;

    double averageNumberOfHits = 0;
    const double averageNumberOfHitsRef = (!isReferenceGeant4Version) ? 3071.17 : 342.37;

    double averageNumberOfTracks = 0;
    const double averageNumberOfTracksRef = (!isReferenceGeant4Version) ? 2161.43 : 10.11;

    TVector3 averagePosition = {};
    const TVector3 averagePositionRef = (!isReferenceGeant4Version) ? TVector3(-38.8987, 27.5536, 91.3969)
                                                                    : TVector3(-17.8046, -32.5019, -31.8353);

    const double tolerance = 0.001;

    for (size_t i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);

        averageTotalEnergy += event->GetTotalDepositedEnergy() / nEvents;
        averageSensitiveEnergy += event->GetSensitiveVolumeEnergy() / nEvents;
        averageNumberOfHits += event->GetNumberOfHits() / nEvents;
        averageNumberOfTracks += event->GetNumberOfTracks() / nEvents;
        averagePosition += event->GetMeanPositionInVolume(0) * (1.0 / nEvents);
    }

    cout << "Average total energy: " << averageTotalEnergy << " keV" << endl;
    cout << "Average sensitive energy: " << averageSensitiveEnergy << " keV" << endl;
    cout << "Average number of hits: " << averageNumberOfHits << endl;
    cout << "Average number of tracks: " << averageNumberOfTracks << endl;
    cout << "Average position: (" << averagePosition.x() << ", " << averagePosition.y() << ", "
         << averagePosition.z() << ") mm" << endl;

    if (TMath::Abs(averageNumberOfHits - averageNumberOfHitsRef) / averageNumberOfHitsRef > tolerance) {
        cout << "The average number of hits does not match the reference value of " << averageNumberOfHitsRef
             << endl;
        return 8;
    }

    if (TMath::Abs(averageNumberOfTracks - averageNumberOfTracksRef) / averageNumberOfTracksRef > tolerance) {
        cout << "The average number of tracks does not match the reference value of "
             << averageNumberOfTracksRef << endl;
        return 9;
    }

    if (TMath::Abs(averageSensitiveEnergy - averageSensitiveEnergyRef) / averageSensitiveEnergyRef >
        tolerance) {
        cout << "The average sensitive volume energy does not match the reference value of "
             << averageSensitiveEnergyRef << endl;
        return 10;
    }

    if (TMath::Abs(averageTotalEnergy - averageTotalEnergyRef) / averageTotalEnergyRef > tolerance) {
        cout << "The average total energy does not match the reference value of " << averageTotalEnergyRef
             << endl;
        return 11;
    }

    if (TMath::Abs(averagePosition.Mag() - averagePositionRef.Mag()) / averagePositionRef.Mag() > tolerance) {
        cout << "The average position does not match the reference value of "
             << "(" << averagePositionRef.x() << ", " << averagePositionRef.y() << ", "
             << averagePositionRef.z() << ") mm" << endl;
        return 12;
    }

    return 0;
}