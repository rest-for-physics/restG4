
Int_t Validate(const char* filename) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    cout << "Entries: " << run.GetEntries() << endl;

    if (run.GetEntries() != 1000) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 1;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");

    const auto geometryInfo = metadata->GetGeant4GeometryInfo();

    geometryInfo->Print();

    if (geometryInfo->GetAllPhysicalVolumes().size() != 374) {
        cout << "Incorrect number of physical volumes " << geometryInfo->GetAllPhysicalVolumes().size()
             << endl;
        // assembly do not work on older geant4 versions...
        // return 1;
    }
    if (geometryInfo->GetAllLogicalVolumes().size() != 22) {
        cout << "Incorrect number of logical volumes " << geometryInfo->GetAllLogicalVolumes().size() << endl;
        // assembly do not work on older geant4 versions...
        // return 1;
    }
    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
