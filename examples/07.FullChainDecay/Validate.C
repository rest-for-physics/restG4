
Int_t Validate(string filename, Int_t nDaughters) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    auto event = (TRestGeant4Event*)run.GetInputEvent();

    std::vector<string> eventTags;
    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);
        eventTags.emplace_back((string)event->GetSubEventTag());
    }

    std::sort(eventTags.begin(), eventTags.end());

    auto iter = std::unique(eventTags.begin(), eventTags.end());

    eventTags.erase(iter, eventTags.end());

    cout << "Daughter isotopes: " << eventTags.size() << endl;
    for (const auto& tag : eventTags) {
        cout << tag << " ";
    }
    cout << endl;

    if (eventTags.size() != nDaughters) {
        cout << "Wrong number of isotopes found! " << eventTags.size() << " found vs " << nDaughters
             << " required" << endl;
        return 13;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
