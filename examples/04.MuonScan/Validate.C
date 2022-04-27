
Int_t Validate(string fname) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun* run = new TRestRun(fname);

    if (run->GetParentRunNumber() != 0) {
        cout << "Parent run number value : " << run->GetParentRunNumber() << endl;
        cout << "The parent run number from restG4 generated file should be 0" << endl;
        return 1;
    }

    if (run->GetRunNumber() != 111) {
        cout << "Run number value : " << run->GetRunNumber() << endl;
        cout << "The run number on the Muon example should be 111 by default!" << endl;
        return 2;
    }

    if (run->GetRunType() != "restG4") {
        cout << "Run type : " << run->GetRunType() << endl;
        cout << "The run type of restG4 generated data should be 'restG4'!" << endl;
        return 3;
    }

    if (run->GetRunTag() != "MuonsFromPoint") {
        cout << "Run tag : " << run->GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'MuonsFromPoint'!" << endl;
        return 4;
    }

    if (run->GetEntries() != 10) {
        cout << "Run entries : " << run->GetEntries() << endl;
        cout << "The Muon simulation should always generate 10 events." << endl;
        return 5;
    }

    cout << "Testing reading of Geant4 metadata class" << endl;
    TRestGeant4Metadata* g4Md = (TRestGeant4Metadata*)run->GetMetadataClass("TRestGeant4Metadata");
    if (!g4Md) {
        cout << "Problem reading Geant4 metadata class!" << endl;
        return 6;
    }
    g4Md->PrintMetadata();

    if (g4Md->GetNumberOfActiveVolumes() != 2) {
        cout << "The number of registered volumes is not 2!" << endl;
        return 7;
    }

    TRestGeant4Event* ev = (TRestGeant4Event*)run->GetInputEvent();
    run->GetEntry(9);

    cout << "Total energy : " << ev->GetTotalDepositedEnergy() << endl;
    Int_t en = (Int_t)(100 * ev->GetTotalDepositedEnergy());
    cout << "Energy integer : " << en << endl;
    cout << "Sensitive volume energy : " << ev->GetSensitiveVolumeEnergy() << endl;
    cout << "Number of hits : " << ev->GetNumberOfHits() << endl;
    cout << "Number of tracks : " << ev->GetNumberOfTracks() << endl;

    Int_t X = (Int_t)(100 * ev->GetMeanPositionInVolume(0).X());
    Int_t Y = (Int_t)(100 * ev->GetMeanPositionInVolume(0).Y());
    Int_t Z = (Int_t)(100 * ev->GetMeanPositionInVolume(0).Z());

    cout << "x: " << X << " y: " << Y << " z: " << Z << endl;

    if (en != 1206) {
        cout << "Error in total energy" << endl;
        return 8;
    }

    en = (Int_t)(100 * ev->GetSensitiveVolumeEnergy());
    if (en != 496) {
        cout << "Error in total energy" << endl;
        return 9;
    }

    if (ev->GetNumberOfHits() != 102) {
        cout << "Error in the number of hits" << endl;
        return 10;
    }

    if (ev->GetNumberOfTracks() != 1) {
        cout << "Error in the number of tracks" << endl;
        return 11;
    }
    if (X != 0 || Y != 0 || Z != -30180) {
        cout << "Error in the event mean position" << endl;
        return 12;
    }

    cout << "Number of hits in active volume 0: " << ev->GetNumberOfHits(0) << endl;
    cout << "Number of hits in active volume 1: " << ev->GetNumberOfHits(1) << endl;
    if (ev->GetNumberOfHits(0) != 51 || ev->GetNumberOfHits(1) != 70) {
        cout << "Error in number of hits at the active volumes" << endl;
        return 13;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    // Other tests like opening other metadata classes. Detector TGeoManager, etc.

    return 0;
}
