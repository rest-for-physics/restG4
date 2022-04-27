
Int_t Validate() {
    TFile* f = new TFile("plots.root");

    TH1F* h1 = (TH1F*)f->Get("1MeV_5um");
    TH1F* h2 = (TH1F*)f->Get("5MeV_5um");
    TH1F* h3 = (TH1F*)f->Get("5MeV_1um");

    cout << "Entries: h1: " << h1->Integral() << ", h2: " << h2->Integral() << ", h3: " << h3->Integral()
         << endl;


    TRestRun* run = new TRestRun("data/Run_g4Analysis_1MeV_5um.root");
    TRestAnalysisTree* aT = run->GetAnalysisTree();

    run->GetEntry(100);

    Int_t thetaInt = (int)(1000 * aT->GetObservableValue<Double_t>("g4Ana_thetaPrimary"));
    Int_t phiInt = (int)(1000 * aT->GetObservableValue<Double_t>("g4Ana_phiPrimary"));
    Double_t theta = aT->GetObservableAverage("g4Ana_thetaPrimary");
    Double_t phi = aT->GetObservableAverage("g4Ana_phiPrimary");

    cout << "entry 100, thetaInt: " << thetaInt << ", phiInt: " << phiInt << endl;
    cout << "average theta: " << theta << ", average phi: " << phi << endl;

    
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

    if (thetaInt != 515) {
        cout << "Wrong theta angle value for entry 100!" << endl;
        cout << "Theta value is : " << thetaInt << endl;
        return 60;
    }


    if (phiInt != 1539) {
        cout << "Wrong phi angle value for entry 100!" << endl;
        cout << "Phi value is : " << phiInt << endl;
        return 60;
    }


    if ((int)(1000 * theta) != 1589) {
        cout << "Wrong theta angle average!" << endl;
        cout << "Theta angle average : " << theta << " while it should be : " << 1.589 << endl;
        return 80;
    }


    if ((int)(1000 * phi) != 27) {
        cout << "Wrong phi angle average!" << endl;
        cout << "Phi angle average : " << phi << " while it should be : " << 0.027 << endl;
        return 90;
    }


    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
