
#include "SimulationManager.h"

#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4Threading.hh>

using namespace std;

SimulationManager::SimulationManager() {
    // Only master thread should create the SimulationManager
    if (!G4Threading::IsMasterThread()) {
        cout << "Only master thread should create the SimulationManager!" << endl;
        exit(1);
    }
}

void SimulationManager::InitializeOutputManager() { fOutputManager = new OutputManager(this); }

SimulationManager::~SimulationManager() {
    delete fRestRun;
    delete fRestGeant4Metadata;
    delete fRestGeant4PhysicsLists;
}

size_t SimulationManager::InsertEvent(std::unique_ptr<TRestGeant4Event>& event) {
    fEventContainerMutex.lock();
    fEventContainer.push(std::move(event));
    auto size = fEventContainer.size();
    fEventContainerMutex.unlock();
    return size;
}

void SimulationManager::WriteEvents() {
    if (G4Threading::IsMultithreadedApplication()) {
        lock_guard<mutex> guard(fEventContainerMutex);
    }

    if (fEventContainer.empty()) {
        return;
    }

    while (!fEventContainer.empty()) {
        fEvent = *fEventContainer.front();

        const auto eventTree = fRestRun->GetEventTree();
        if (eventTree != nullptr) {
            eventTree->Fill();
        }

        const auto analysisTree = fRestRun->GetAnalysisTree();
        if (analysisTree != nullptr) {
            analysisTree->SetEventInfo(&fEvent);
            analysisTree->Fill();
        }

        fEventContainer.pop();
    }
}

// OutputManager
OutputManager::OutputManager(const SimulationManager* simulationManager)
    : fSimulationManager(const_cast<SimulationManager*>(simulationManager)) {
    // this class should only exist on the threads performing the simulation
    if (G4Threading::IsMasterThread() && G4Threading::IsMultithreadedApplication()) {
        G4cout << "Error in 'OutputManager', this instance should never exist" << endl;
        exit(1);
    }
}

void OutputManager::UpdateEvent() {
    auto event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    fEvent = make_unique<TRestGeant4Event>(event, *fSimulationManager->fRestGeant4Metadata);
}

bool OutputManager::IsEmptyEvent() const { return !fEvent || fEvent->fTrack.empty(); }

bool OutputManager::IsValidEvent() const {
    if (IsEmptyEvent()) return false;
    if (fSimulationManager->fRestGeant4Metadata->GetSaveAllEvents()) return true;
    if (fEvent->GetSensitiveVolumeEnergy() <= 0) return false;
    return true;
}

void OutputManager::FinishAndSubmitEvent() {
    if (IsEmptyEvent()) return;

    if (IsValidEvent()) {
        size_t numberOfInsertedEvents = fSimulationManager->InsertEvent(fEvent);
    }
    UpdateEvent();
}

// Geant4Lib

TRestGeant4Event::TRestGeant4Event(const G4Event* event, const TRestGeant4Metadata& metadata)
    : TRestGeant4Event() {
    SetID(event->GetEventID());
    SetOK(true);
    time_t system_time = time(nullptr);

    SetTime((Double_t)system_time);

    // Defining if the hits in a given volume will be stored
    for (int i = 0; i < metadata.GetNumberOfActiveVolumes(); i++) {
        if (metadata.GetStorageChance(i) >= 1.00) {
            ActivateVolumeForStorage(i);
        } else {
            Double_t randomNumber = G4UniformRand();
            if (metadata.GetStorageChance(i) >= randomNumber) {
                ActivateVolumeForStorage(i);
            } else {
                DisableVolumeForStorage(i);
            }
        }
    }
}

TRestGeant4Track::TRestGeant4Track(const G4Track* track) : TRestGeant4Track() {
    fTrackID = track->GetTrackID();
    fParentID = track->GetParentID();

    auto particle = track->GetParticleDefinition();
    fParticleName = particle->GetParticleName();
    /*
    fParticleID = particle->GetPDGEncoding();
    fParticleType = particle->GetParticleType();
    fParticleSubType = particle->GetParticleSubType();
    */
    if (track->GetCreatorProcess()) {
        fCreatorProcess = track->GetCreatorProcess()->GetProcessName();
    } else {
        fCreatorProcess = "IsPrimaryParticle";
    }

    fKineticEnergy = track->GetKineticEnergy() / CLHEP::keV;

    fWeight = track->GetWeight();

    fGlobalTimestamp = track->GetGlobalTime() / CLHEP::second;

    const G4ThreeVector& trackOrigin = track->GetPosition();
    fTrackOrigin = {trackOrigin.x(), trackOrigin.y(), trackOrigin.z()};
}