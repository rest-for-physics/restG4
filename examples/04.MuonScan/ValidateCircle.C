
Int_t ValidateCircle(string fname) {
    gSystem->Load("/usr/local/rest-for-physics/lib/libRestFramework.so");
    gSystem->Load("/usr/local/rest-for-physics/lib/libRestGeant4.so");

    TRestRun run(fname);
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    if (run.GetRunTag() != "MuonsFromCircle") {
        cout << "Run tag : " << run.GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'MuonsFromCircle'" << endl;
        return 4;
    }

    Double_t rMean = 0;
    Double_t rMin = TMath::Infinity();
    Double_t rMax = 0;
    for (Int_t n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);
        Double_t x = event->GetPrimaryEventOrigin().X();
        Double_t z = event->GetPrimaryEventOrigin().Z();

        const auto r = TMath::Sqrt(x * x + z * z) / 100;  // cm

        rMean += r / run.GetEntries();
        if (r < rMin) {
            rMin = r;
        }
        if (r > rMax) {
            rMax = r;
        }
    }

    if (rMean < 2.75 || rMean > 3.25) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "R_mean (cm): " << rMean << endl;
        return 5;
    }

    if (rMin > 0.5) {
        cout << "The minimum radius of the distribution is wrong!" << endl;
        cout << "R_min (cm): " << rMin << endl;
        return 6;
    }

    if (rMax > 4.0 || rMax < 3.75) {
        cout << "The maximum radius of the distribution is wrong!" << endl;
        cout << "R_max (cm): " << rMax << endl;
        return 7;
    }

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() < 600 || run.GetEntries() > 750) {
        cout << "The number of entries is wrong!" << endl;
        cout << "Number of entries : " << run.GetEntries() << endl;
        return 8;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";

    return 0;
}
