
#include "SimulationManager.h"

#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include <G4Threading.hh>
#include <Randomize.hh>

#include "SteppingAction.h"

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

void SimulationManager::InitializeUserDistributions() {
    auto random = []() { return (double)G4UniformRand(); };

    for (int i = 0; i < fRestGeant4Metadata->GetNumberOfSources(); i++) {
        fRestGeant4Metadata->GetParticleSource(i)->SetRandomMethod(random);
    }

    TRestGeant4ParticleSource* source = fRestGeant4Metadata->GetParticleSource(0);

    if (TRestGeant4PrimaryGeneratorTypes::StringToEnergyDistributionTypes(
            source->GetEnergyDistributionType().Data()) ==
        TRestGeant4PrimaryGeneratorTypes::EnergyDistributionTypes::TH1D) {
        TFile file(source->GetEnergyDistributionFilename());
        auto distribution = (TH1D*)file.Get(source->GetEnergyDistributionNameInFile());

        if (!distribution) {
            RESTError << "Error when trying to find energy spectrum" << RESTendl;
            RESTError << "File: " << source->GetEnergyDistributionFilename() << RESTendl;
            RESTError << "Spectrum name: " << source->GetEnergyDistributionNameInFile() << RESTendl;
            exit(1);
        }

        fPrimaryEnergyDistribution = *distribution;
    }

    if (TRestGeant4PrimaryGeneratorTypes::StringToAngularDistributionTypes(
            source->GetAngularDistributionType().Data()) ==
        TRestGeant4PrimaryGeneratorTypes::AngularDistributionTypes::TH1D) {
        TFile file(source->GetAngularDistributionFilename());
        auto distribution = (TH1D*)file.Get(source->GetAngularDistributionNameInFile());

        if (!distribution) {
            RESTError << "Error when trying to find angular spectrum" << RESTendl;
            RESTError << "File: " << source->GetAngularDistributionFilename() << RESTendl;
            RESTError << "Spectrum name: " << source->GetAngularDistributionNameInFile() << RESTendl;
            exit(1);
        }

        fPrimaryAngularDistribution = *distribution;
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

    // initialize active volume lookup set
    const auto metadata = fSimulationManager->GetRestMetadata();
    for (size_t i = 0; i < metadata->GetNumberOfActiveVolumes(); i++) {
        const TString& activeVolume = metadata->GetActiveVolumeName(i);
        fActiveVolumes.insert(activeVolume.Data());
    }
}

void OutputManager::UpdateEvent() {
    auto event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    fEvent = make_unique<TRestGeant4Event>(event);
    fEvent->InitializeReferences(fSimulationManager->GetRestRun());
}

bool OutputManager::IsEmptyEvent() const { return !fEvent || fEvent->fTracks.empty(); }

bool OutputManager::IsValidEvent() const {
    if (IsEmptyEvent()) return false;
    if (fSimulationManager->GetRestMetadata()->GetSaveAllEvents()) return true;
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
        const auto& lastTrack = fEvent->fTracks.back();
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

void OutputManager::RecordStep(const G4Step* step) { fEvent->InsertStep(step); }

void OutputManager::AddSensitiveEnergy(Double_t energy, const char* physicalVolumeName) {
    fEvent->AddEnergyToSensitiveVolume(energy);
    /*
        const TString physicalVolumeNameNew = fSimulationManager->GetRestMetadata()->GetGeant4GeometryInfo()
                                                  .GetAlternativeNameFromGeant4PhysicalName(physicalVolumeName);
                                                  */
}

void OutputManager::AddEnergyToVolumeForProcess(Double_t energy, const char* volumeName,
                                                const char* processName) {
    if (energy <= 0) {
        return;
    }
    if (fEvent->fEnergyInVolumePerProcess[volumeName].count(processName) == 0) {
        fEvent->fEnergyInVolumePerProcess[volumeName][processName] = 0;
    }
    fEvent->fEnergyInVolumePerProcess[volumeName][processName] += energy;

    fEvent->fTotalDepositedEnergy += energy;
}

// Geant4Lib

TRestGeant4Event::TRestGeant4Event(const G4Event* event) : TRestGeant4Event() {
    SetID(event->GetEventID());
    SetOK(true);
    time_t system_time = time(nullptr);

    SetTime((Double_t)system_time);

    auto primaryVertex = event->GetPrimaryVertex();
    const auto& position = primaryVertex->GetPosition();
    fPrimaryPosition = {position.x() / CLHEP::mm, position.y() / CLHEP::mm, position.z() / CLHEP::mm};
    for (int i = 0; i < primaryVertex->GetNumberOfParticle(); i++) {
        const auto& primaryParticle = primaryVertex->GetPrimary(i);
        fPrimaryParticleNames.emplace_back(primaryParticle->GetParticleDefinition()->GetParticleName());
        fPrimaryEnergies.emplace_back(primaryParticle->GetKineticEnergy() / CLHEP::keV);
        const auto& momentum = primaryParticle->GetMomentumDirection();
        fPrimaryDirections.emplace_back(momentum.x(), momentum.y(), momentum.z());
    }

    return;

    // TODO: move this
    // Defining if the hits in a given volume will be stored
    const auto metadata = GetGeant4Metadata();
    for (int i = 0; i < metadata->GetNumberOfActiveVolumes(); i++) {
        if (metadata->GetStorageChance(i) >= 1.00) {
            // ActivateVolumeForStorage(i);
        } else {
            Double_t randomNumber = G4UniformRand();
            if (metadata->GetStorageChance(i) >= randomNumber) {
                // ActivateVolumeForStorage(i);
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

    if (fTracks.empty() && IsSubEvent()) {
        // First track of sub-event (primary)
        fSubEventPrimaryParticleName = track->GetParticleDefinition()->GetParticleName();
        fSubEventPrimaryEnergy = track->GetKineticEnergy() / CLHEP::keV;
        const auto& position = track->GetPosition();
        fSubEventPrimaryPosition = {position.x() / CLHEP::mm, position.y() / CLHEP::mm,
                                    position.z() / CLHEP::mm};
        const auto& momentum = track->GetMomentumDirection();
        fSubEventPrimaryDirection = {momentum.x(), momentum.y(), momentum.z()};
    }

    fTrackIDToTrackIndex[track->GetTrackID()] = fTracks.size();  // before insertion

    fTracks.emplace_back(track);

    auto& insertedTrack = fTracks.back();

    insertedTrack.SetHits(fInitialStep);
    insertedTrack.SetEvent(this);

    TRestGeant4Track* parentTrack = GetTrackByID(track->GetParentID());
    if (parentTrack) {
        parentTrack->AddSecondaryTrackID(track->GetTrackID());
    }

    return true;
}

void TRestGeant4Event::UpdateTrack(const G4Track* track) { fTracks.back().UpdateTrack(track); }

void TRestGeant4Event::InsertStep(const G4Step* step) {
    if (step->GetTrack()->GetCurrentStepNumber() == 0) {
        // initial step (from SteppingVerbose) is generated before TrackingAction can insert the first track
        fInitialStep = TRestGeant4Hits();
        fInitialStep.SetEvent(this);
        fInitialStep.InsertStep(step);
    } else {
        fTracks.back().InsertStep(step);
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
    if (track->GetCreatorProcess() != nullptr) {
        fCreatorProcess = track->GetCreatorProcess()->GetProcessName();
    } else {
        fCreatorProcess = "PrimaryGenerator";
    }

    fInitialKineticEnergy = track->GetKineticEnergy() / CLHEP::keV;

    fWeight = track->GetWeight();

    fGlobalTimestamp = track->GetGlobalTime() / CLHEP::second;

    const G4ThreeVector& trackOrigin = track->GetPosition();
    fInitialPosition = {trackOrigin.x(), trackOrigin.y(), trackOrigin.z()};
}

void TRestGeant4Track::InsertStep(const G4Step* step) { fHits.InsertStep(step); }

void TRestGeant4Track::UpdateTrack(const G4Track* track) {
    if (track->GetTrackID() != fTrackID) {
        G4cout << "Geant4Track::UpdateTrack - mismatch of trackID!" << endl;
        exit(1);
    }

    fLength = track->GetTrackLength() / CLHEP::mm;
    fTimeLength = track->GetGlobalTime() / CLHEP::second - fGlobalTimestamp;
}

Int_t TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(const G4VProcess* process) {
    return process->GetProcessType() * 1000 + process->GetProcessSubType();
}

void TRestGeant4Hits::InsertStep(const G4Step* step) {
    const G4Track* track = step->GetTrack();

    TRestGeant4Metadata* metadata = GetGeant4Metadata();

    const auto& geometryInfo = metadata->GetGeant4GeometryInfo();

    // Variables that describe a step are taken.
    const auto& volumeName = geometryInfo.GetAlternativeNameFromGeant4PhysicalName(
        (TString &&) step->GetPreStepPoint()->GetPhysicalVolume()->GetName());

    if (!SimulationManager::GetOutputManager()->IsActiveVolume(volumeName) &&
        step->GetTrack()->GetCurrentStepNumber() != 0) {
        // we always store the first step
        return;
    }

    const auto& particle = step->GetTrack()->GetDefinition();
    const auto& particleID = particle->GetPDGEncoding();
    const auto& particleName = particle->GetParticleName();

    metadata->fGeant4PhysicsInfo.InsertParticleName(particleID, particleName);

    const auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4String processName = "Init";
    G4String processTypeName = "Init";
    Int_t processID = 0;
    if (track->GetCurrentStepNumber() != 0) {
        // 0 = Init step (G4SteppingVerbose) process is not defined for this step
        processName = process->GetProcessName();
        processTypeName = G4VProcess::GetProcessTypeName(process->GetProcessType());
        processID = TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(process);
    }

    metadata->fGeant4PhysicsInfo.InsertProcessName(processID, processName);

    const auto energy = step->GetTotalEnergyDeposit() / CLHEP::keV;
    const auto trackKineticEnergy = step->GetTrack()->GetKineticEnergy() / CLHEP::keV;

    auto sensitiveVolumeName =
        geometryInfo.GetAlternativeNameFromGeant4PhysicalName(metadata->GetSensitiveVolume());

    if (particle->GetParticleName() == "geantino" && sensitiveVolumeName.Data() == volumeName) {
        metadata->SetSaveAllEvents(true);
    }

    G4Track* aTrack = step->GetTrack();

    Double_t x = aTrack->GetPosition().x() / CLHEP::mm;
    Double_t y = aTrack->GetPosition().y() / CLHEP::mm;
    Double_t z = aTrack->GetPosition().z() / CLHEP::mm;

    const TVector3 hitPosition(x, y, z);
    const Double_t hitGlobalTime = step->GetPreStepPoint()->GetGlobalTime() / CLHEP::second;
    const G4ThreeVector& momentum = step->GetPreStepPoint()->GetMomentumDirection();

    AddHit(hitPosition, energy, hitGlobalTime);  // this increases fNHits

    fProcessID.emplace_back(processID);
    fVolumeID.emplace_back(geometryInfo.GetIDFromVolume(volumeName));
    fKineticEnergy.emplace_back(step->GetPreStepPoint()->GetKineticEnergy() / CLHEP::keV);
    fMomentumDirection.emplace_back(momentum.x(), momentum.y(), momentum.z());

    SimulationManager::GetOutputManager()->AddEnergyToVolumeForProcess(energy, volumeName, processName);
}
