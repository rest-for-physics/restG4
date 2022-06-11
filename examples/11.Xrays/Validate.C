#include "TRestRun.h"

Int_t Validate(const char* inputFile) {
    RestRun run(inputFile);
    TRestAnalysisTree* aTree = run.GetAnalysisTree();

    TH1D* eDepHist = new TH1D("eDepHist", "Energy deposited in detector", 150, 0, 15);
    TH1D* primaryHist = new TH1D("primaryHist", "Energy of generated primary", 150, 0, 15);

    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);

        Double_t eDep = aTree->GetObservableValue<Double_t>("g4Ana_totalEDep");
        Double_t ePrimary = aTree->GetObservableValue<Double_t>("g4Ana_energyPrimary");

        if (eDep > 0) eDepHist->Fill(eDep);
        primaryHist->Fill(ePrimary);
    }

    FILE* f = fopen("output.txt", "wt");
    for (int n = 0; n < eDepHist->GetNbinsX(); n++) {
        Double_t edepCounts = eDepHist->GetBinContent(n + 1);
        Double_t primaryCounts = primaryHist->GetBinContent(n + 1);

        Double_t eff = edepCounts / primaryCounts;
        Double_t err = TMath::Sqrt(edepCounts) / primaryCounts;
        fprintf(f, "%4.2f\t%5.0lf\t%5.0lf\t%4.4lf\t%4.4lf\n", eDepHist->GetBinCenter(n + 1), edepCounts,
                primaryCounts, eff, err);
    }
    fclose(f);

    // Plot
    TCanvas* c1 = new TCanvas("c1", "response matrix", 1200, 800);
    c1->GetFrame()->SetBorderSize(6);
    c1->GetFrame()->SetBorderMode(-1);

    TPad* pad1 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
    pad1->Divide(2, 2);
    pad1->Draw();

    TH2D* matrix_hist = new TH2D("matrix hist", "response matrix", 150, 0, 15, 150, 0, 15);
    matrix_hist->GetXaxis()->SetTitle("deposited Energy [keV]");
    matrix_hist->GetYaxis()->SetTitle("primary Energy [keV]");
    gPad->SetLogz();
    // gPad->SetTheta(90);
    // gPad->SetPhi(0);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.15);

    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);

        Double_t eDep = aTree->GetObservableValue<Double_t>("g4Ana_totalEDep");
        Double_t ePrimary = aTree->GetObservableValue<Double_t>("g4Ana_energyPrimary");

        if (eDep > 0) matrix_hist->Fill(eDep, ePrimary);
    }

    TH1D* prHist = matrix_hist->ProjectionY("PrimaryEnergy");
    pad1->cd(1);
    prHist->Draw();

    TH1D* depHist = matrix_hist->ProjectionX("DepositEnergy");
    pad1->cd(2);
    depHist->Draw();

    pad1->cd(3);
    primaryHist->Draw();

    pad1->cd(4);
    eDepHist->Draw();

    //   gPad->cd(3);
    //   matrix_hist->Draw("SURF1");

    c1->Print("montecarlo.png");

    TH1D* effHisto = new TH1D("Efficiency", "Argon isobutane at 1 bar", 150, 0, 15);
    for (int n = 1; n <= effHisto->GetNbinsX(); n++)
        effHisto->SetBinContent(n, eDepHist->GetBinContent(n) / primaryHist->GetBinContent(n));

    TCanvas* c2 = new TCanvas("c2", "Efficiency", 800, 600);
    c2->GetFrame()->SetBorderSize(6);
    c2->GetFrame()->SetBorderMode(-1);

    effHisto->SetStats(0);
    effHisto->Draw();

    c2->Print("efficiency.png");

    cout << "Integral: " << eDepHist->Integral() << endl;
    if (eDepHist->Integral() < 500) {
        cout << "Something went wrong. Number of counts inside deposits integral too low" << endl;
        return 1;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
