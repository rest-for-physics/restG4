
Int_t Validate(string fname) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun* run = new TRestRun(fname);

    TRestGeant4Event* ev = (TRestGeant4Event*)run->GetInputEvent();
    run->GetEntry(0);

    Int_t nHitsEv0 = ev->GetNumberOfHits();
    cout << "Nhits : " << nHitsEv0 << endl;
    if (nHitsEv0 < 20 || nHitsEv0 > 40) {
        cout << "The number of hits in the recoil is not in the 20-40 range" << endl;
        return 14;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
