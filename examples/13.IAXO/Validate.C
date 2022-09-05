
Int_t Validate(const char* filename) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    cout << "Entries: " << run.GetEntries() << endl;

    if (run.GetEntries() != 1) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 1;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");

    const auto geometryInfo = metadata->GetGeant4GeometryInfo();

    geometryInfo.Print();

    if (geometryInfo.GetAllPhysicalVolumes().size() != 374) {
        cout << "Incorrect number of physical volumes " << geometryInfo.GetAllPhysicalVolumes().size()
             << endl;
    }
    if (geometryInfo.GetAllLogicalVolumes().size() != 22) {
        cout << "Incorrect number of logical volumes " << geometryInfo.GetAllLogicalVolumes().size() << endl;
    }

    // event information
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    run.GetEntry(0);

    if (event->GetNumberOfTracks() != 28059) {
        cout << "Incorrect number of tracks: " << event->GetNumberOfTracks() << endl;
    }
    if (event->GetNumberOfHits() != 87189) {
        cout << "Incorrect number of hits: " << event->GetNumberOfHits() << endl;
    }

    if (TMath::Abs(event->GetSensitiveVolumeEnergy() - 72.436883) > 1e-4) {
        cout << "Incorrect sensitive volume energy: " << event->GetSensitiveVolumeEnergy() << endl;
    }

    const auto scintillatorEnergy = event->GetEnergyInVolume(
        "VetoSystem_vetoSystemBack_vetoLayerBack1_assembly-17.veto3_scintillatorVolume-800.0mm-f1a5df6c");
    if (TMath::Abs(scintillatorEnergy - 1348.3876) > 1e-4) {
        cout << "Incorrect scintillator volume energy: " << scintillatorEnergy << endl;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";

    return 0;
}
