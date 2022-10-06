Int_t ValidateResponse(std::string fName) {
    TRestRun run(fName);
    TRestAnalysisTree* aTree = run.GetAnalysisTree();
    aTree->Draw("g4Ana_totalEdep");
    TH1D* h = (TH1D*)aTree->GetHistogram();

    TRestGeant4Metadata* md = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");

    Double_t efficiency = h->Integral() / md->GetNumberOfEvents();

    Double_t efficiencyReference = 0.91;
    std::cout << "Overall efficiency : " << efficiency << std::endl;

    if (TMath::Abs(efficiency - efficiencyReference) > 0.03) {
        cout << "The efficiency does not match the reference value of " << efficiencyReference << endl;
        return 1;
    }

    return 0;
}
