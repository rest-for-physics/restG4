
Int_t Validate() {
    TFile* f = new TFile("plots.root");

    TH1F* h1 = (TH1F*)f->Get("1MeV_5um");
    if (h1->Integral() != 1578) {
        cout << "Wrong number of alphas produced at 1MeV_5um!" << endl;
        cout << "Histogram contains : " << h1->Integral() << endl;
        return 10;
    }

    TH1F* h2 = (TH1F*)f->Get("5MeV_5um");
    if (h2->Integral() != 7182) {
        cout << "Wrong number of alphas produced at 5MeV_5um!" << endl;
        cout << "Histogram contains : " << h2->Integral() << endl;
        return 20;
    }

    TH1F* h3 = (TH1F*)f->Get("5MeV_1um");
    if (h3->Integral() != 9421) {
        cout << "Wrong number of alphas produced at 5MeV_1um!" << endl;
        cout << "Histogram contains : " << h3->Integral() << endl;
        return 30;
    }

    TRestRun* run = new TRestRun("data/Run_g4Analysis_1MeV_5um.root");
    TRestAnalysisTree* aT = run->GetAnalysisTree();

    // Testing an arbitrary entry
    run->GetEntry(100);

    Int_t thetaInt = (int)(1000 * aT->GetObservableValue<Double_t>("g4Ana_thetaPrimary"));
    if (thetaInt != 515) {
        cout << "Wrong theta angle value for entry 100!" << endl;
        return 60;
    }

    Int_t phiInt = (int)(1000 * aT->GetObservableValue<Double_t>("g4Ana_phiPrimary"));
    if (phiInt != 1539) {
        cout << "Wrong theta angle value for entry 100!" << endl;
        return 60;
    }

    Double_t theta = aT->GetObservableAverage("g4Ana_thetaPrimary");
    if ((int)(1000 * theta) != 1589) {
        cout << "Wrong theta angle average!" << endl;
        cout << "Theta angle average : " << theta << " while it should be : " << 1.589 << endl;
        return 80;
    }

    Double_t phi = aT->GetObservableAverage("g4Ana_phiPrimary");
    if ((int)(1000 * phi) != 27) {
        cout << "Wrong phi angle average!" << endl;
        cout << "Phi angle average : " << phi << " while it should be : " << 0.027 << endl;
        return 90;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
