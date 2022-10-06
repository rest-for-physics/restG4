Int_t ValidateFe55(std::string fName) {
    TRestRun run(fName);

    TRestAnalysisTree* aTree = run.GetAnalysisTree();

    aTree->Draw("g4Ana_totalEdep");

    TH1D* h = (TH1D*)aTree->GetHistogram();

    Int_t peakCounts = h->Integral(h->FindFixBin(5.7), h->FindFixBin(5.9));

    if (peakCounts < 1600 || peakCounts > 1750) {
        std::cout << "Problem on Fe55 gamma peak identification." << std::endl;
        std::cout << "Peak counts found : " << peakCounts << std::endl;
        std::cout << "The result should be in the range ( 1600, 1750 )" << std::endl;
    }

    return 0;
}
