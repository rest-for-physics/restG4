
Int_t Validate(string fname) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun* run = new TRestRun(fname);

    TRestAnalysisTree* aT = run->GetAnalysisTree();
    aT->Draw("g4Ana_totalEdep>>h(100,1,10)");

    TH1D* histo = (TH1D*)gDirectory->Get("h");

    Double_t mean = histo->GetMean();
    cout << "Mean : " << histo->GetMean() << endl;

    if (mean > 8.5 || mean < 7.5) {
        cout << "Error copper fluorescences are not observed" << endl;
        return 13;
    }

    Int_t entries = histo->GetEntries();
    cout << "Entries : " << entries << endl;

    if (entries < 8000) {
        cout << "The number of gammas observed is too low" << endl;
        return 14;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
