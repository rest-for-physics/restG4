//
// Created by lobis on 10/11/2021.
//

#include "GlobalManager.h"

#include <TRestGDMLParser.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestRun.h>
#include <TTree.h>

#include <G4RunManager.hh>
#include <G4Threading.hh>

#include "spdlog/spdlog.h"

using namespace std;

GlobalManager* GlobalManager::pinstance_ = nullptr;

GlobalManager* GlobalManager::Instance() {
    if (!pinstance_) {
        pinstance_ = new GlobalManager();
    }
    return pinstance_;
}

GlobalManager::GlobalManager()
    : fRestGeant4Metadata(nullptr),
      fRestRun(nullptr),
      fRestGeant4PhysicsLists(nullptr),
      fRestGDMLParser(new TRestGDMLParser) {
    // Master thread should create the GlobalManager, worker threads should spawn after the manager is created
    // by the master thread
    if (!G4Threading::IsMasterThread()) {
        cout
            << "GlobalManager::GlobalManager - ERROR GlobalManager should never be created by a worker thread"
            << endl;
        exit(1);
    }
}

GlobalManager::~GlobalManager() {
    cout << "GlobalManager::~GlobalManager (Destructor)" << endl;

    delete fRestGeant4Metadata;
    delete fRestGeant4PhysicsLists;
    delete fRestRun;

    //
    delete fPrimaryEnergyDistribution;
    delete fPrimaryAngularDistribution;
}

void GlobalManager::InitializeFromConfigFile(const TString& rmlFile) {
    if (fRestGeant4Metadata || fRestRun || fRestGeant4PhysicsLists) {
        cout << "GlobalManager::InitializeRestGeant4Metadata - ERROR rest classes should not be initialized "
                "twice"
             << endl;
        exit(1);
    }

    fInputConfigFile = rmlFile;

    InitializeRestGeant4Metadata(fInputConfigFile);
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(fRestGeant4Metadata->GetSeed());

    InitializeRestGeant4PhysicsLists(fInputConfigFile);

    InitializeRestRun(fInputConfigFile);

    InitializeTrees();
}

void GlobalManager::InitializeRestGeant4Metadata(const TString& rmlFile) {
    fRestGeant4Metadata = new TRestGeant4Metadata(const_cast<char*>(rmlFile.Data()));

    REST_Verbose_Level verboseLevel = fRestGeant4Metadata->GetVerboseLevel();
    if (verboseLevel == REST_Verbose_Level::REST_Debug) {
        spdlog::set_level(spdlog::level::debug);
    } else if (verboseLevel == REST_Verbose_Level::REST_Info) {
        spdlog::set_level(spdlog::level::info);
    } else if (verboseLevel == REST_Verbose_Level::REST_Essential) {
        spdlog::set_level(spdlog::level::warn);
    } else if (verboseLevel == REST_Verbose_Level::REST_Silent) {
        spdlog::set_level(spdlog::level::err);
    }

    spdlog::set_pattern("[%T][%^%l%$][thread %t]: %v");

    string geant4Version = TRestTools::Execute("geant4-config --version");
    fRestGeant4Metadata->SetGeant4Version(geant4Version);
    // GDML geometry parsing
    // This call will generate a new single file GDML output
    fRestGDMLParser->Load((string)fRestGeant4Metadata->Get_GDML_Filename());

    // We redefine the value of the GDML file to be used in DetectorConstructor.
    fRestGeant4Metadata->Set_GDML_Filename(fRestGDMLParser->GetOutputGDMLFile());
    fRestGeant4Metadata->SetGeometryPath("");

    fRestGeant4Metadata->Set_GDML_Reference(fRestGDMLParser->GetGDMLVersion());
    fRestGeant4Metadata->SetMaterialsReference(fRestGDMLParser->GetEntityVersion("materials"));

    SetSaveAllEventsFlag(fRestGeant4Metadata->GetSaveAllEvents());

    // Primary
    if (fRestGeant4Metadata->GetParticleSource(0)->GetEnergyDistType() == "TH1D") {
        TString fileFullPath = (TString)fRestGeant4Metadata->GetParticleSource(0)->GetSpectrumFilename();
        TFile file(fileFullPath);
        TString spectrumName = fRestGeant4Metadata->GetParticleSource(0)->GetSpectrumName();

        auto primaryEnergyDistribution = (TH1D*)file.Get(spectrumName);
        fPrimaryEnergyDistribution =
            (TH1D*)
                primaryEnergyDistribution->Clone();  // This MUST be static cast! not dynamic! (not sure why)

        if (!fPrimaryEnergyDistribution) {
            spdlog::error(
                "GlobalManager::InitializeRestGeant4Metadata - Error initializing energy distribution "
                "spectrum ({}) from file {}",
                spectrumName, fileFullPath);
            exit(1);
        } else {
            fPrimaryEnergyDistribution->SetDirectory(nullptr);

            cout << "ENERGY DISTRIBUTION: " << fPrimaryEnergyDistribution->GetName() << " "
                 << fPrimaryEnergyDistribution << endl;

            fPrimaryEnergyDistributionMin = fRestGeant4Metadata->GetParticleSource(0)->GetMinEnergy();
            if (fPrimaryEnergyDistributionMin < 0) {
                fPrimaryEnergyDistributionMin = 0;
            }
            fPrimaryEnergyDistributionMax = fRestGeant4Metadata->GetParticleSource(0)->GetMaxEnergy();
            if (fPrimaryEnergyDistributionMax < 0) {
                fPrimaryEnergyDistributionMax = 0;
            }
        }
    }

    if (fRestGeant4Metadata->GetParticleSource(0)->GetAngularDistType() == "TH1D") {
        TString fileFullPath = (TString)fRestGeant4Metadata->GetParticleSource(0)->GetAngularFilename();

        TFile file(fileFullPath);

        TString spectrumName = fRestGeant4Metadata->GetParticleSource(0)->GetAngularName();

        auto primaryAngularDistribution = (TH1D*)file.Get(spectrumName);
        fPrimaryAngularDistribution =
            (TH1D*)
                primaryAngularDistribution->Clone();  // This MUST be static cast! not dynamic! (not sure why)

        if (!fPrimaryAngularDistribution) {
            spdlog::error(
                "GlobalManager::InitializeRestGeant4Metadata - Error initializing angular distribution "
                "spectrum ({}) from file {}",
                spectrumName, fileFullPath);
            exit(1);
        } else {
            fPrimaryAngularDistribution->SetDirectory(nullptr);

            cout << "ANGULAR DISTRIBUTION: " << fPrimaryAngularDistribution->GetName() << " "
                 << fPrimaryAngularDistribution << endl;
        }
    }
}

void GlobalManager::InitializeRestRun(const TString& rmlFile) {
    fRestRun = new TRestRun();
    fRestRun->LoadConfigFromFile(const_cast<char*>(rmlFile.Data()));
    TString runTag = fRestRun->GetRunTag();
    if (runTag == "Null" || runTag == "") {
        fRestRun->SetRunTag(fRestGeant4Metadata->GetTitle());
    }

    fRestRun->SetRunType("restG4");

    fRestRun->AddMetadata(fRestGeant4Metadata);
    fRestRun->AddMetadata(fRestGeant4PhysicsLists);

    fRestRun->PrintMetadata();

    fRestRun->FormOutputFile();
}

void GlobalManager::InitializeRestGeant4PhysicsLists(const TString& rmlFile) {
    fRestGeant4PhysicsLists = new TRestGeant4PhysicsLists(const_cast<char*>(rmlFile.Data()));
}

size_t GlobalManager::InsertEvent(std::unique_ptr<TRestGeant4Event>& event) {
    fEventContainerMutex.lock();
    fEventContainer.push(std::move(event));
    auto size = fEventContainer.size();
    fEventContainerMutex.unlock();
    return size;
}

void GlobalManager::InitializeTrees() {
    // Event Tree
    fEventTree = fRestRun->GetEventTree();
    // fEventTree->Branch("EventBranch", &fEvent);
    fRestRun->AddEventBranch(&fEvent);

    fAnalysisTree = fRestRun->GetAnalysisTree();
}

Long64_t GlobalManager::GetEntries() {
    if (!fEventTree) return 0;
    if (G4Threading::IsMultithreadedApplication()) {
        lock_guard<mutex> guard(fEventContainerMutex);
    }
    return fEventTree->GetEntries();
}

void GlobalManager::FillEvents() {
    /*
    spdlog::debug("GlobalManager::WriteEvents");

    if (!fSaveFlag) {
        spdlog::debug("GlobalManager::FillEvents() - Saving events is disabled");
        return;
    }
    */
    if (G4Threading::IsMultithreadedApplication()) {
        lock_guard<mutex> guard(fEventContainerMutex);
    }

    if (fEventContainer.empty()) {
        return;
    }

    while (!fEventContainer.empty()) {
        fEvent = *fEventContainer.front();
        fEventTree->Fill();
        fAnalysisTree->Fill();
        fEventContainer.pop();
    }

    // fEventTree->Write();
    /*
    const auto beforeNumberOfEvents = fEventTree->GetEntries();

    spdlog::debug("GlobalManager::FillEvents - Saved {} events into {} (total {} events)",
                  fEventTree->GetEntries() - beforeNumberOfEvents, fFile->GetName(),
                  fEventTree->GetEntries());
                  */
}

void GlobalManager::WriteEvents() { fEventTree->Write(); }
