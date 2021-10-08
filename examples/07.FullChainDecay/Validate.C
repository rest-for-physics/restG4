
Int_t Validate(string fname, Int_t nDaughters) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun* run = new TRestRun(fname);

    TRestGeant4Event* ev = (TRestGeant4Event*)run->GetInputEvent();

    std::vector<string> evTag;
    for (int n = 0; n < run->GetEntries(); n++) {
        run->GetEntry(n);
        evTag.push_back((string)ev->GetSubEventTag());
    }

    std::sort(evTag.begin(), evTag.end());

    auto iter = std::unique(evTag.begin(), evTag.end());

    evTag.erase(iter, evTag.end());

    cout << "Daughter isotopes: " << evTag.size() << endl;
    for (unsigned int n = 0; n < evTag.size(); n++) cout << evTag[n] << " ";
    cout << endl;

    if (evTag.size() != nDaughters) {
        cout << "Wrong number of isotopes found!" << endl;
        return 13;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
