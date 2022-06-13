
Int_t Validate(const char* filename, Int_t nDaughters) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    std::set<string> eventTagsUnique;
    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);
        eventTagsUnique.insert((string)event->GetSubEventTag());
    }

    cout << "Daughter isotopes: " << eventTagsUnique.size() << endl;
    for (const auto& tag : eventTagsUnique) {
        cout << tag << " ";
    }
    cout << endl;

    if (eventTagsUnique.size() != nDaughters) {
        cout << "Wrong number of isotopes found! " << eventTagsUnique.size() << " vs. reference value of "
             << nDaughters << endl;
        return 1;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
