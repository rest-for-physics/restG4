
Int_t Validate() {
    TFile* file = TFile::Open("plots.root");

    TH1F* h1 = file->Get<TH1F>("1MeV_5um");
    TH1F* h2 = file->Get<TH1F>("5MeV_5um");
    TH1F* h3 = file->Get<TH1F>("5MeV_1um");

    cout << "Entries: h1: " << h1->Integral() << ", h2: " << h2->Integral() << ", h3: " << h3->Integral()
         << endl;

    TRestRun run("data/Run_g4Analysis_1MeV_5um.root");
    TRestAnalysisTree* analysisTree = run.GetAnalysisTree();

    run.GetEntry(100);

    Double_t thetaSample = analysisTree->GetObservableValue<Double_t>("g4Ana_thetaPrimary");
    constexpr Double_t thetaSampleRef = 2.85884;

    Double_t phiSample = analysisTree->GetObservableValue<Double_t>("g4Ana_phiPrimary");
    constexpr Double_t phiSampleRef = -1.59582;

    Double_t thetaAverage = analysisTree->GetObservableAverage("g4Ana_thetaPrimary");
    constexpr Double_t thetaAverageRef = 1.5767;

    Double_t phiAverage = analysisTree->GetObservableAverage("g4Ana_phiPrimary");
    constexpr Double_t phiAverageRef = 0.0719113;

    cout << "entry 100, theta: " << thetaSample << ", phi: " << phiSample << endl;
    cout << "average theta: " << thetaAverage << ", average phi: " << phiAverage << endl;

    if (h1->Integral() != 1926) {
        cout << "Wrong number of alphas produced at (h1) 1MeV_5um!" << endl;
        cout << "Histogram contains : " << h1->Integral() << endl;
        return 10;
    }

    if (h2->Integral() != 7487) {
        cout << "Wrong number of alphas produced at (h2) 5MeV_5um!" << endl;
        cout << "Histogram contains : " << h2->Integral() << endl;
        return 20;
    }

    if (h3->Integral() != 9432) {
        cout << "Wrong number of alphas produced at (h3) 5MeV_1um!" << endl;
        cout << "Histogram contains : " << h3->Integral() << endl;
        return 30;
    }

    if (TMath::Abs(thetaSample - thetaSampleRef) / thetaSampleRef >= 0.01) {
        cout << "Wrong theta angle value for entry 100!" << endl;
        cout << "Theta value is : " << thetaSample << " while it should be : " << thetaSampleRef << endl;
        return 60;
    }

    if (TMath::Abs(phiSample - phiSampleRef) / phiSampleRef >= 0.01) {
        cout << "Wrong phi angle value for entry 100!" << endl;
        cout << "Phi value is : " << phiSample << " while it should be : " << phiSampleRef << endl;
        return 60;
    }

    if (TMath::Abs(thetaAverage - thetaAverageRef) / thetaAverageRef >= 0.01) {
        cout << "Wrong theta angle average!" << endl;
        cout << "Theta angle average : " << thetaAverage << " while it should be : " << thetaAverageRef
             << endl;
        return 80;
    }

    if (TMath::Abs(phiAverage - phiAverageRef) / phiAverageRef >= 0.01) {
        cout << "Wrong phi angle average!" << endl;
        cout << "Phi angle average : " << phiAverage << " while it should be : " << phiAverageRef << endl;
        return 90;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
