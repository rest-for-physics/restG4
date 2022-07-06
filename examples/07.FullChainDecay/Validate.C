
Int_t Validate(const char* filename, Int_t nDaughters) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    std::set<string> uniqueSubEventPrimaryParticleNames;
    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);
        if (!event->IsSubEvent()) {
            continue;
        }
        uniqueSubEventPrimaryParticleNames.insert(event->GetSubEventPrimaryEventParticleName().Data());
    }

    cout << "Daughter isotopes: " << uniqueSubEventPrimaryParticleNames.size() << endl;
    for (const auto& primary : uniqueSubEventPrimaryParticleNames) {
        cout << primary << " ";
    }
    cout << endl;

    if (uniqueSubEventPrimaryParticleNames.size() != nDaughters) {
        cout << "Wrong number of isotopes found! " << uniqueSubEventPrimaryParticleNames.size()
             << " vs. reference value of " << nDaughters << endl;
        return 1;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
