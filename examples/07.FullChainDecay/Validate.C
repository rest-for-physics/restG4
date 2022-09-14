
Int_t Validate(const char* filename, Int_t nDaughters) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    std::set<string> uniqueDecayingParticleNames;
    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);
        if (!event->IsSubEvent()) {
            uniqueDecayingParticleNames.insert(event->GetPrimaryEventParticleName().Data());
        } else {
            uniqueDecayingParticleNames.insert(event->GetSubEventPrimaryEventParticleName().Data());
        }
    }

    cout << "Decaying isotopes: " << uniqueDecayingParticleNames.size() << endl;
    for (const auto& primary : uniqueDecayingParticleNames) {
        cout << primary << " ";
    }
    cout << endl;

    if (uniqueDecayingParticleNames.size() != nDaughters) {
        cout << "Wrong number of isotopes found! " << uniqueDecayingParticleNames.size()
             << " vs. reference value of " << nDaughters << endl;
        return 1;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
