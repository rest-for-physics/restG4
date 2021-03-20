
// A basic example that extracts the Mean Free Path of the recoils travelling in the medium
TH1D* RecoilMFP(string fname) {
    TRestRun* run = new TRestRun(fname);
    TRestGeant4Event* ev = (TRestGeant4Event*)run->GetInputEvent();

    Int_t nEntries = run->GetEntries();
    cout << "Number of entries: " << nEntries << endl;

    // Creating a histogram with 1000 bins in 1mm. I.e. 1um step.
    TH1D* hist = new TH1D("Recoil MFP", "Recoil MFP", 1000, 0, 1);

    for (Int_t n = 0; n < nEntries; n++) {
        run->GetEntry(n);

        Int_t mTracks = ev->GetNumberOfTracks();
        for (Int_t m = 0; m < mTracks; m++) {
            TString pName = ev->GetTrack(m)->GetParticleName();
            Double_t Ekin = ev->GetTrack(m)->GetKineticEnergy();

            /// Track origin (If using the recoil example they are launched anyway from (0,0,0)
            Double_t xO = ev->GetTrack(m)->GetTrackOrigin().X();
            Double_t yO = ev->GetTrack(m)->GetTrackOrigin().Y();
            Double_t zO = ev->GetTrack(m)->GetTrackOrigin().Z();

            /// Last hit position
            Int_t lHits = ev->GetTrack(m)->GetNumberOfHits();
            Double_t xR = ev->GetTrack(m)->GetHits()->GetPosition(lHits - 1).X();
            Double_t yR = ev->GetTrack(m)->GetHits()->GetPosition(lHits - 1).Y();
            Double_t zR = ev->GetTrack(m)->GetHits()->GetPosition(lHits - 1).Z();

            TVector3 vector(xO - xR, yO - yR, zO - zR);

            cout << "Particle : " << pName << " Kinetic energy : " << Ekin << endl;
            vector.Print();

            // We take only those tracks that fulfill that are a F20 recoil. Other conditions?
            // Condition in the energy to generate a histo for a given particle range?
            if (pName == "F20") {
                // We get the recoil calculated distance and fill the histogram
                hist->Fill(vector.Mag());
            }
        }

        // Now we create a canvas and Draw the histogram?
        // TCanvas c;
        // hist->Draw();
    }
    return hist;
}
