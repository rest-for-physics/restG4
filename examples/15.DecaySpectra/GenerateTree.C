#include <filesystem>
#include <memory>

using namespace std;

void GenerateTree(const char* filenameOrDirectoryWithFiles) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    vector<string> filenames;
    if (filesystem::is_directory(filenameOrDirectoryWithFiles)) {
        for (const auto& entry : filesystem::directory_iterator(filenameOrDirectoryWithFiles)) {
            filenames.push_back(entry.path());
            cout << "Adding file: " << entry.path() << endl;
        }
    } else {
        filenames.push_back(filenameOrDirectoryWithFiles);
    }

    TFile* outputFile = nullptr;
    const map<string, string> particleToString = {
        {"gamma", "gammas"},
        /*
        {"e-",      "electrons"},
        {"e+",      "positrons"},
        {"proton",  "protons"},
        {"neutron", "neutrons"},
        {"alpha",   "alphas"},
         */
    };
    vector<string> particlesToRecord;
    for (const auto& particle : particleToString) {
        particlesToRecord.push_back(particle.first);
    }
    map<string, TTree*> trees;
    map<string, double> energies;

    bool first = true;
    for (const auto& filename : filenames) {
        cout << "Processing file: " << filename << endl;

        TRestRun run(filename);
        cout << "Number of entries: " << run.GetEntries() << endl;

        TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
        run.GetEntry(0);

        if (first) {
            const TString primaryName = event->GetPrimaryEventParticleName();
            // Create a new tree to store gammas from the primary particle
            outputFile =
                new TFile(TString::Format("ParticlesFrom%sDecay.root", primaryName.Data()), "RECREATE");

            // create a new tree for each particle
            for (const auto& particle : particlesToRecord) {
                // we cannot use the particle name (e.g. gamma) as a tree name, because it is a reserved word
                trees[particle] =
                    new TTree(particleToString.at(particle).c_str(), particleToString.at(particle).c_str());
            }

            // create a branch for each particle
            for (const auto& particle : particlesToRecord) {
                energies[particle] = 0;
                trees[particle]->Branch("energy", &energies[particle], "energy/D");
            }
            first = false;
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
                        // We can use either the kinetic energy at the point of exit from the detector or the
                        // initial energy depending on our goals energies[particleName] =
                        // hits.GetKineticEnergy(i);
                        energies[particleName] = track.GetInitialKineticEnergy();
                        trees[particleName]->Fill();
                    }
                }
            }
        }
    }

    outputFile->Write();
    outputFile->Close();

    /*
     * Histograms can be merged into a single file using the following root code:

TFile* U238 = TFile::Open("U238/gammas.root");
TFile* Th232 = TFile::Open("Th232/gammas.root");

TH1D* U238h = U238->Get<TH1D>("histogram");
TH1D* Th232h = Th232->Get<TH1D>("histogram");

// rename histograms
U238h->SetName("U238");
Th232h->SetName("Th232");

// write to file
TFile* outputFile = new TFile("merged.root", "RECREATE");

U238h->Write();
Th232h->Write();

outputFile->Close();

     */
}
