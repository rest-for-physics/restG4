#include "TRestRun.h"

Int_t Validate(const char* inputFile) {
    TRestRun run(inputFile);
    TRestAnalysisTree* aTree = run.GetAnalysisTree();

    auto eDepHist = new TH1D("eDepHist", "Energy deposited in detector", 150, 0, 15);
    auto primaryHist = new TH1D("primaryHist", "Energy of generated primary", 150, 0, 15);

    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);

        Double_t eDep = aTree->GetObservableValue<Double_t>("g4Ana_totalEDep");
        Double_t ePrimary = aTree->GetObservableValue<Double_t>("g4Ana_energyPrimary");

        if (eDep > 0) eDepHist->Fill(eDep);
        primaryHist->Fill(ePrimary);
    }

    FILE* f = fopen("output.txt", "wt");
    for (int n = 0; n < eDepHist->GetNbinsX(); n++) {
        Double_t eDepCounts = eDepHist->GetBinContent(n + 1);
        Double_t primaryCounts = primaryHist->GetBinContent(n + 1);

        Double_t eff = eDepCounts / primaryCounts;
        Double_t err = TMath::Sqrt(eDepCounts) / primaryCounts;
        fprintf(f, "%4.2f\t%5.0lf\t%5.0lf\t%4.4lf\t%4.4lf\n", eDepHist->GetBinCenter(n + 1), eDepCounts,
                primaryCounts, eff, err);
    }
    fclose(f);

    // Plot
    auto c1 = new TCanvas("c1", "response matrix", 1200, 800);
    c1->GetFrame()->SetBorderSize(6);
    c1->GetFrame()->SetBorderMode(-1);

    auto pad = new TPad("pad", "This is pad", 0.01, 0.02, 0.99, 0.97);
    pad->Divide(2, 2);
    pad->Draw();

    auto matrixHist = new TH2D("matrix hist", "response matrix", 150, 0, 15, 150, 0, 15);
    matrixHist->GetXaxis()->SetTitle("deposited Energy [keV]");
    matrixHist->GetYaxis()->SetTitle("primary Energy [keV]");
    gPad->SetLogz();
    // gPad->SetTheta(90);
    // gPad->SetPhi(0);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.15);

    for (int n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);

        Double_t eDep = aTree->GetObservableValue<Double_t>("g4Ana_totalEDep");
        Double_t ePrimary = aTree->GetObservableValue<Double_t>("g4Ana_energyPrimary");

        if (eDep > 0) matrixHist->Fill(eDep, ePrimary);
    }

    TH1D* prHist = matrixHist->ProjectionY("PrimaryEnergy");
    pad->cd(1);
    prHist->Draw();

    TH1D* depHist = matrixHist->ProjectionX("DepositEnergy");
    pad->cd(2);
    depHist->Draw();

    pad->cd(3);
    primaryHist->Draw();

    pad->cd(4);
    eDepHist->Draw();

    // gPad->cd(3);
    // matrixHist->Draw("SURF1");

    c1->Print("montecarlo.png");

    auto effHisto = new TH1D("Efficiency", "Argon isobutane at 1 bar", 150, 0, 15);
    for (int n = 1; n <= effHisto->GetNbinsX(); n++)
        effHisto->SetBinContent(n, eDepHist->GetBinContent(n) / primaryHist->GetBinContent(n));

    auto c2 = new TCanvas("c2", "Efficiency", 800, 600);
    c2->GetFrame()->SetBorderSize(6);
    c2->GetFrame()->SetBorderMode(-1);

    effHisto->SetStats(0);
    effHisto->Draw();

    c2->Print("efficiency.png");

    constexpr double eDepIntegralReference = 500.0;
    cout << "Integral: " << eDepHist->Integral() << endl;
    if (eDepHist->Integral() < eDepIntegralReference) {
        cout << "Number of counts inside deposits integral " << eDepHist->Integral()
             << " does not match reference value of " << eDepIntegralReference << endl;
        return 1;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
