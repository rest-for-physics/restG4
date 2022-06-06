
Int_t Validate(string fname) {
    TRestRun* run = new TRestRun(fname);

    if (run->GetParentRunNumber() != 0) {
        cout << "Parent run number value : " << run->GetParentRunNumber() << endl;
        cout << "The parent run number from restG4 generated file should be 0" << endl;
        return 1;
    }

    if (run->GetRunNumber() != 1) {
        cout << "Run number value : " << run->GetRunNumber() << endl;
        cout << "The run number on the validation chain should be 1 by default!" << endl;
        return 2;
    }

    if (run->GetRunType() != "restG4") {
        cout << "Run type : " << run->GetRunType() << endl;
        cout << "The run type of restG4 generated data should be 'restG4'!" << endl;
        return 3;
    }

    if (run->GetRunTag() != "NLDBD") {
        cout << "Run tag : " << run->GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'NLDBD'!" << endl;
        return 4;
    }

    if (run->GetEntries() != 100) {
        cout << "Run entries : " << run->GetEntries() << endl;
        cout << "The NLDBD simulation is launched from gas. It should always generate 100 events." << endl;
        return 5;
    }

    cout << "Testing reading of Geant4 metadata class" << endl;
    TRestGeant4Metadata* geant4Metadata = (TRestGeant4Metadata*)run->GetMetadataClass("TRestGeant4Metadata");
    if (!geant4Metadata) {
        cout << "Problem reading Geant4 metadata class!" << endl;
        return 6;
    }
    geant4Metadata->PrintMetadata();

    if (geant4Metadata->GetNumberOfActiveVolumes() != 1) {
        cout << "The number of registered volumes is not 1!" << endl;
        return 7;
    }

    TRestGeant4Event* event = run->GetInputEvent<TRestGeant4Event>();
    run->GetEntry(9);

    cout << "Total energy : " << event->GetTotalDepositedEnergy() << endl;
    Int_t en = (Int_t)(100 * event->GetTotalDepositedEnergy());
    cout << "Energy integer : " << en << endl;
    cout << "Sensitive volume energy : " << event->GetSensitiveVolumeEnergy() << endl;
    cout << "Number of hits : " << event->GetNumberOfHits() << endl;
    cout << "Number of tracks : " << event->GetNumberOfTracks() << endl;

    Int_t X = (Int_t)(100 * event->GetMeanPositionInVolume(0).X());
    Int_t Y = (Int_t)(100 * event->GetMeanPositionInVolume(0).Y());
    Int_t Z = (Int_t)(100 * event->GetMeanPositionInVolume(0).Z());

    cout << "x: " << X << " y: " << Y << " z: " << Z << endl;

    if (en != 245700) {
        cout << "Error in total energy" << endl;
        return 8;
    }

    if (event->GetNumberOfHits() != 401) {
        cout << "Error in the number of hits" << endl;
        return 10;
    }

    if (event->GetNumberOfTracks() != 11) {
        cout << "Error in the number of tracks" << endl;
        return 11;
    }
    if (X != -52709 || Y != 37592 || Z != 59585) {
        cout << "Error in the event mean position" << endl;
        return 12;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    // Other tests like opening other metadata classes. Detector TGeoManager, etc.

    return 0;
}
