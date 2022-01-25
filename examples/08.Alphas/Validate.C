
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

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
