
Int_t Validate(string filename) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun* run = new TRestRun(filename);

    cout << "Entries: " << run->GetEntries() << endl;

    if (run->GetEntries() != 1000) {
        cout << "Bad number of entries: " << run->GetEntries() << endl;
        return -1;
    }

    auto metadata = (TRestGeant4Metadata*)run->GetMetadataClass("TRestGeant4Metadata");

    const auto geometryInfo = metadata->GetGeometryInfo();

    geometryInfo->Print();

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
