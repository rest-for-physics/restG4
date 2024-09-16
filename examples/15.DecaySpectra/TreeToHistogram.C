
using namespace std;

void TreeToHistogram(const char* filename) {
    TFile* file = TFile::Open(filename);

    TTree* tree = nullptr;
    tree = file->Get<TTree>("gammas");

    TH1D* histogram = new TH1D("histogram", "histogram", 10000, 0, 5000);

    double energy;
    tree->SetBranchAddress("energy", &energy);
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        histogram->Fill(energy);
    }
    // write to file
    TFile* outputFile = new TFile("gammas.root", "RECREATE");

    histogram->Write();
    outputFile->Close();
}
