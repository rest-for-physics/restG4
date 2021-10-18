
Int_t ValidateWall(string fname) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun* run = new TRestRun(fname);
    TRestGeant4Event* ev = (TRestGeant4Event*)run->GetInputEvent();

    if (run->GetRunTag() != "MuonsFromWall") {
        cout << "Run tag : " << run->GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'MuonsFromWall!" << endl;
        return 4;
    }

    Double_t r = 0;
    for (Int_t n = 0; n < run->GetEntries(); n++) {
        run->GetEntry(n);
        Double_t x = ev->GetPrimaryEventOrigin().X();
        Double_t y = ev->GetPrimaryEventOrigin().Y();

        r += x * x + y * y;
    }
    r /= run->GetEntries();

    if ((Int_t)(1000. * r) / 1000 < 10000 || (Int_t)(1000. * r) / 1000 > 20000) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "R: " << (Int_t)(1000. * r) / 1000 << endl;
        return 5;
    }

    cout << "Run entries: " << run->GetEntries() << endl;
    if (run->GetEntries() < 350 || run->GetEntries() > 450) {
        cout << "The number of entries is not between 350 and 450!" << endl;
        cout << "Number of entries : " << run->GetEntries() << endl;
        return 6;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    // Other tests like opening other metadata classes. Detector TGeoManager, etc.

    return 0;
}
