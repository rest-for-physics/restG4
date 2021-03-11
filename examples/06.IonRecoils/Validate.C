
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

    Double_t duration = run->GetEndTimestamp() - run->GetStartTimestamp();

    cout << "Duration : " << duration << endl;
    if (duration <= 5 || duration > 200) {
        cout << "The duration of the run in seconds is invalid" << endl;
        return 15;
    }

    cout.precision(12);

    cout << "Time since simulation ended : " << time(NULL) - run->GetEndTimeStamp() << endl;
    Double_t delay = time(NULL) - run->GetEndTimeStamp();

    if (delay < 0 || delay > 20) {
        cout << "The end timestamp is probably wrong" << endl;
        return 16;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
