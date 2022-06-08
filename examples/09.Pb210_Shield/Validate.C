
Int_t Validate(const char* filename) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    constexpr double tolerance = 0.005;
    constexpr int numberOfEntriesRef = 252;

    if (TMath::Abs(run.GetEntries() - numberOfEntriesRef) / numberOfEntriesRef > tolerance) {
        cout << "Number of entries: " << run.GetEntries() << "does not match the reference value of "
             << numberOfEntriesRef << endl;
        return 1;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
