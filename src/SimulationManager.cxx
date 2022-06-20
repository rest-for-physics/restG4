
#include "SimulationManager.h"

#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include <G4Threading.hh>

using namespace std;

thread_local OutputManager* SimulationManager::fOutputManager = nullptr;

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

void OutputManager::RecordTrack(const G4Track* track) {
    if (!IsValidTrack(track)) {
        return;
    }
    fEvent->InsertTrack(track);

    if (fEvent->fSubEventID > 0) {
        const auto& lastTrack = fEvent->fTrack.back();
        assert(lastTrack.GetTrackID() == track->GetTrackID());
        // TODO
        /*
        bool isSubEventPrimary = fEvent->IsTrackSubEventPrimary(lastTrack.fTrackID);
        if (isSubEventPrimary) {
            spdlog::debug(
                "OutputManager::RecordTrack - Setting track ID {} as SubEventPrimaryTrack of EventID {} "
                "(SubEventID {}). Track info: {} - Created by "
                "{} - ParentID: {}",
                lastTrack.fTrackID, fEvent->fEventID, fEvent->fSubEventID, lastTrack.fParticleName,
                lastTrack.fCreatorProcess, lastTrack.fParentID);
        }
         */
    }
}

void OutputManager::UpdateTrack(const G4Track* track) {
    if (!IsValidTrack(track)) {
        return;
    }
    fEvent->UpdateTrack(track);
}

void OutputManager::RecordStep(const G4Step* step) {
    fEvent->InsertStep(step, *fSimulationManager->fRestGeant4Metadata);
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

bool TRestGeant4Event::InsertTrack(const G4Track* track) {
    if (fInitialStep.GetNumberOfHits() != 1) {
        cout << "fInitialStep does not have exactly one step! Problem with stepping verbose" << endl;
        exit(1);
    }
    fTrack.emplace_back(track);
    fTrack.back().SetHits(fInitialStep);

    fTrackIDToTrackIndex[track->GetTrackID()] = fTrack.size() - 1;

    if (fTrack.empty()) {
        /*
        // primary for the sub-event
        fSubEventPrimaryParticleName = track->GetParticleDefinition()->GetParticleName();
        fSubEventPrimaryEnergy = track->GetKineticEnergy() / CLHEP::keV;
        const auto& position = track->GetPosition();
        fSubEventPrimaryPosition =
            TVector3(position.x() / CLHEP::mm, position.y() / CLHEP::mm, position.z() / CLHEP::mm);
        const auto& momentum = track->GetMomentumDirection();
        fSubEventPrimaryMomentum =
            TVector3(momentum.x() / CLHEP::mm, momentum.y() / CLHEP::mm, momentum.z() / CLHEP::mm);
            */
    }
    return true;
}

void TRestGeant4Event::UpdateTrack(const G4Track* track) { fTrack.back().UpdateTrack(track); }

void TRestGeant4Event::InsertStep(const G4Step* step, TRestGeant4Metadata& metadata) {
    if (step->GetTrack()->GetCurrentStepNumber() == 0) {
        // initial step (from SteppingVerbose) is generated before TrackingAction can insert the first track
        fInitialStep = TRestGeant4Hits();
        fInitialStep.InsertStep(step, metadata);
    } else {
        fTrack.back().InsertStep(step, metadata);
    }
}

bool OutputManager::IsValidTrack(const G4Track*) const { return true; }

bool OutputManager::IsValidStep(const G4Step*) const { return true; }

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

void TRestGeant4Track::InsertStep(const G4Step* step, TRestGeant4Metadata& metadata) {
    fHits.InsertStep(step, metadata);
}

void TRestGeant4Track::UpdateTrack(const G4Track* track) {
    if (track->GetTrackID() != fTrackID) {
        G4cout << "Geant4Track::UpdateTrack - mismatch of trackID!" << endl;
        exit(1);
    }

    fTrackLength = track->GetTrackLength() / CLHEP::mm;

    // auto steppingAction = (SteppingAction*)G4EventManager::GetEventManager()->GetUserSteppingAction();
    // auto secondaries = steppingAction->GetfSecondary();
    // fNumberOfSecondaries = (int)secondaries->size();
}

void TRestGeant4Hits::InsertStep(const G4Step* step, TRestGeant4Metadata& metadata) {
    const G4Track* track = step->GetTrack();

    const auto& geometryInfo = metadata.GetGeant4GeometryInfo();
    // Variables that describe a step are taken.
    const auto& volumeName = geometryInfo.GetAlternativeNameFromGeant4PhysicalName(
        (TString &&) step->GetPreStepPoint()->GetPhysicalVolume()->GetName());

    const auto& particle = step->GetTrack()->GetDefinition();
    const auto& particleID = particle->GetPDGEncoding();
    const auto& particleName = particle->GetParticleName();

    metadata.fGeant4PhysicsInfo.InsertParticleName(particleID, particleName);

    const auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
    const auto& processID = process->GetProcessType() * 1000 + process->GetProcessSubType();
    const auto& processName = process->GetProcessName();

    metadata.fGeant4PhysicsInfo.InsertProcessName(processID, processName);

    const auto totalEnergyDeposit = step->GetTotalEnergyDeposit() / CLHEP::keV;
    const auto trackKineticEnergy = step->GetTrack()->GetKineticEnergy() / CLHEP::keV;

    auto sensitiveVolumeName =
        geometryInfo.GetAlternativeNameFromGeant4PhysicalName(metadata.GetSensitiveVolume());

    if (particle->GetParticleName() == "geantino" && sensitiveVolumeName.Data() == volumeName) {
        metadata.SetSaveAllEvents(true);
    }

    if (!step->GetPostStepPoint()->GetProcessDefinedStep()) {
        G4cout << endl;
        G4cout << endl;
        G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        G4cout << "An ERROR was detected on the G4Step process definition." << endl;
        G4cout << "This is a sign of problem in the restG4 particle definition" << endl;
        G4cout << endl;
        G4cout << "E.g. A definition of a gamma with 0keV energy" << endl;
        G4cout << endl;
        G4cout << "Please, review your TRestGeant4Metadata RML definition" << endl;
        G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        G4cout << endl;
        G4cout << endl;

        exit(0);
    }

    G4Track* aTrack = step->GetTrack();

    Double_t x = aTrack->GetPosition().x() / CLHEP::mm;
    Double_t y = aTrack->GetPosition().y() / CLHEP::mm;
    Double_t z = aTrack->GetPosition().z() / CLHEP::mm;

    if (metadata.GetSensitiveVolume() == volumeName) {
        // restG4Event->AddEnergyToSensitiveVolume(totalEnergyDeposit);
    }

    const TVector3 hitPosition(x, y, z);
    const Double_t hitGlobalTime = step->GetPreStepPoint()->GetGlobalTime() / CLHEP::second;
    const G4ThreeVector& momentum = step->GetPreStepPoint()->GetMomentumDirection();
    const TVector3 momentumDirection = TVector3(momentum.x(), momentum.y(), momentum.z());  //.Unit();
}
