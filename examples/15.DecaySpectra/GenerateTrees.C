
using namespace std;

void GenerateTrees(const char* filename) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);
    cout << "Number of entries: " << run.GetEntries() << endl;

    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    run.GetEntry(0);

    const TString primaryName = event->GetPrimaryEventParticleName();
    // Create a new tree to store gammas from the primary particle
    TFile* f = new TFile(TString::Format("ParticlesFrom%sDecay.root", primaryName.Data()), "RECREATE");

    // store particles and name of the root TTree
    const map<string, string> particleToString = {
        {"gamma", "photons"},  {"e-", "electrons"},     {"e+", "positrons"},
        {"proton", "protons"}, {"neutron", "neutrons"}, {"alpha", "alphas"},
    };
    vector<string> particlesToRecord;
    for (const auto& particle : particleToString) {
        particlesToRecord.push_back(particle.first);
    }
    // create a new tree for each particle
    map<string, TTree*> trees;
    for (const auto& particle : particlesToRecord) {
        // we cannot use the particle name (e.g. gamma) as a tree name, because it is a reserved word
        cout << "Creating tree for " << particle << " with name " << particleToString.at(particle) << endl;
        trees[particle] =
            new TTree(particleToString.at(particle).c_str(), particleToString.at(particle).c_str());
    }

    // create a branch for each particle
    map<string, double> energies;
    for (const auto& particle : particlesToRecord) {
        cout << "Creating branch for " << particle << endl;
        energies[particle] = 0;
        trees[particle]->Branch("energy", &energies[particle], "energy/D");
    }

    const string detectorName = "detector";
    cout << "Detector name: " << detectorName << endl;
    cout << "Starting loop over events..." << endl;
    for (int n = 0; n < run.GetEntries(); n++) {
        // print the progress every 1%
        if (n % (run.GetEntries() / 100) == 0) {
            cout << "Progress: " << n * 100 / run.GetEntries() << "%" << endl;
        }
        run.GetEntry(n);
        // iterate over all particles in the event
        for (const auto& track : event->GetTracks()) {
            // if particle not in the list of particles to record, skip it
            const string particleName = track.GetParticleName().Data();
            if (find(particlesToRecord.begin(), particlesToRecord.end(), particleName) ==
                particlesToRecord.end()) {
                continue;
            }
            // record energy at the point of exit from the "detector" sphere (if it exits)
            const auto& hits = track.GetHits();
            for (int i = 0; i < hits.GetNumberOfHits(); i++) {
                if (hits.GetVolumeName(i) == detectorName && hits.GetProcessName(i) == "Transportation") {
                    // this is more accurate than track.GetInitialKineticEnergy()
                    // because it takes into account energy loss in the detector
                    energies[particleName] = hits.GetKineticEnergy(i);
                    trees[particleName]->Fill();
                }
            }
        }
    }

    f->Write();
    f->Close();
}
