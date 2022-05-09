
Int_t ValidateWall(string fname) {
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
        Double_t y = event->GetPrimaryEventOrigin().Y();

        const auto r = TMath::Sqrt(x * x + y * y) / 100;  // cm

        rMean += r / run.GetEntries();
        if (r < rMin) {
            rMin = r;
        }
        if (r > rMax) {
            rMax = r;
        }
    }

    if (rMean < 0.8 || rMean > 1.2) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "R_mean (cm): " << rMean << endl;
        return 5;
    }

    if (rMin > 0.15) {
        cout << "The minimum radius of the distribution is wrong!" << endl;
        cout << "R_min (cm): " << rMin << endl;
        return 6;
    }

    if (rMax > 1.5 || rMax < 1.35) {
        cout << "The maximum radius of the distribution is wrong!" << endl;
        cout << "R_max (cm): " << rMax << endl;
        return 7;
    }

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() < 350 || run.GetEntries() > 450) {
        cout << "The number of entries is not between 350 and 450!" << endl;
        cout << "Number of entries : " << run.GetEntries() << endl;
        return 8;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";

    return 0;
}
