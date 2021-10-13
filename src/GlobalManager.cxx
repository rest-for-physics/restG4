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

    fSaveAllEventsFlag = fRestGeant4Metadata->GetSaveAllEvents();
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

size_t GlobalManager::InsertEvent(std::unique_ptr<TRestGeant4DataEvent>& event) {
    fEventContainerMutex.lock();
    fEventContainer.push(std::move(event));
    auto size = fEventContainer.size();
    fEventContainerMutex.unlock();
    return size;
}

void GlobalManager::InitializeTrees() {
    // Event Tree
    // fEventTree = new TTree("EventTree", "Event Tree");
    fEventTree = fRestRun->GetEventTree();
    fEventTree->Branch("fEvent", &fEvent);

    fAnalysisTree = fRestRun->GetAnalysisTree();
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
