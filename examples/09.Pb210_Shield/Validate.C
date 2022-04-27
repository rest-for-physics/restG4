
Int_t Validate(string fname) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun* run = new TRestRun(fname);

    if (run->GetEntries() != 28) {
        cout << "Entries: " << run->GetEntries() << endl;
        cout << "There was a problem simulation Pb210. The number of entries should be 28" << endl;
        return 10;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
