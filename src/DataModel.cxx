
#include <G4Event.hh>
#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include <G4Threading.hh>
#include <Randomize.hh>

#include "SimulationManager.h"
#include "SteppingAction.h"

static double globalTimeOffset = 0;

using namespace std;

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

    /*
    // TODO: move this
    // Defining if the hits in a given volume will be stored
    const auto metadata = GetGeant4Metadata();
    if (metadata != nullptr) {
        for (int i = 0; i < metadata->GetNumberOfActiveVolumes(); i++) {
            if (metadata->GetStorageChance(i) >= 1.00) {
                ActivateVolumeForStorage(i);
            } else {
                Double_t randomNumber = G4UniformRand();
                if (metadata->GetStorageChance(i) >= randomNumber) {
                    ActivateVolumeForStorage(i);
                } else {
                    DisableVolumeForStorage(i);
                }
            }
        }
    }
    */
}

bool TRestGeant4Event::InsertTrack(const G4Track* track) {
    if (fInitialStep.GetNumberOfHits() != 1) {
        cout << "fInitialStep does not have exactly one step! Problem with stepping verbose" << endl;
        exit(1);
    }

    if ((fTracks.empty() && IsSubEvent()) ||
        (fTracks.empty() && !IsSubEvent() && GetGeant4Metadata()->GetNumberOfSources() == 1)) {
        // First track of sub-event (primary)
        fSubEventPrimaryParticleName = track->GetParticleDefinition()->GetParticleName();
        fSubEventPrimaryEnergy = track->GetKineticEnergy() / CLHEP::keV;
        const auto& position = track->GetPosition();
        fSubEventPrimaryPosition = {position.x() / CLHEP::mm, position.y() / CLHEP::mm,
                                    position.z() / CLHEP::mm};
        const auto& momentum = track->GetMomentumDirection();
        fSubEventPrimaryDirection = {momentum.x(), momentum.y(), momentum.z()};
    }

    if (fTracks.empty() && fSubEventID == 0) globalTimeOffset = 0;

    fTrackIDToTrackIndex[track->GetTrackID()] = int(fTracks.size());  // before insertion

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

    fGlobalTimestamp = track->GetGlobalTime() / CLHEP::microsecond;

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
    fTimeLength = track->GetLocalTime() / CLHEP::microsecond;

    const auto metadata = GetGeant4Metadata();

    if (!metadata->isGlobalTimeReset()) return;

    SetTimeOffset(globalTimeOffset);

    const double precision = metadata->GetResetTimePrecision() / CLHEP::microsecond;
    const double globalTime = track->GetGlobalTime() / CLHEP::microsecond;

    const auto processName = track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    // ResetGloblal time in case of long Radioactive decay
    if (processName == "RadioactiveDecay")
        if (globalTime + precision == globalTime) {
            auto secondaries = track->GetStep()->GetSecondary();
            size_t nSeco = secondaries->size();
            globalTimeOffset += globalTime;
            if (nSeco > 0)
                for (size_t i = 0; i < nSeco; i++) {
                    G4Track* tck = (G4Track*)(*secondaries)[i];
                    tck->SetGlobalTime(0.);  // Reset global time for secondaries
                }
        }
}

Int_t TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(const G4VProcess* process) {
    return process->GetProcessType() * 1000 + process->GetProcessSubType();
}

void TRestGeant4Hits::InsertStep(const G4Step* step) {
    const G4Track* track = step->GetTrack();

    TRestGeant4Metadata* metadata = GetGeant4Metadata();

    const auto& geometryInfo = metadata->GetGeant4GeometryInfo();

    // Get the full name (path) of the physical volume which uniquely identifies it
    auto th = step->GetPreStepPoint()->GetTouchable();
    G4int depth = th->GetHistoryDepth();
    G4String geant4path = "";
    if (depth == 0) {  // it is the world volume
        geant4path = th->GetVolume()->GetName();
    }
    for (G4int i = 1; i <= depth; ++i) {  // start from 1 to skip world volume
        // Move the touchable to level i (0 = current volume, depth = world)
        G4VPhysicalVolume* pv = th->GetVolume(depth - i);
        if (pv) {
            if (geant4path != "") {
                geant4path += geometryInfo.GetPathSeparator().Data();
            }
            geant4path += pv->GetName();
        }
    }
    // convert to the names used in gdml (due to assemblies)
    const auto volumeName = geometryInfo.GetAlternativePathFromGeant4Path(geant4path);

    if (!metadata->IsActiveVolume(volumeName) && step->GetTrack()->GetCurrentStepNumber() != 0) {
        // we always store the first step
        return;
    }

    const bool kill = metadata->IsKillVolume(volumeName);

    const auto& particle = step->GetTrack()->GetDefinition();
    const auto& particleID = particle->GetPDGEncoding();
    const auto& particleName = particle->GetParticleName();

    auto energy = step->GetTotalEnergyDeposit() / CLHEP::keV;

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

    if (kill) {
        processName = "REST-for-physics-kill";
        processTypeName = "REST-for-physics";
        processID = 1000000;  // use id out of range!
        energy = 0;

        step->GetTrack()->SetTrackStatus(fStopAndKill);
    }

    metadata->fGeant4PhysicsInfo.InsertProcessName(processID, processName, processTypeName);

    auto sensitiveVolumeName =
        geometryInfo.GetAlternativeNameFromGeant4PhysicalName(metadata->GetSensitiveVolume());

    G4Track* aTrack = step->GetTrack();

    Double_t x = aTrack->GetPosition().x() / CLHEP::mm;
    Double_t y = aTrack->GetPosition().y() / CLHEP::mm;
    Double_t z = aTrack->GetPosition().z() / CLHEP::mm;

    const TVector3 hitPosition(x, y, z);
    const Double_t hitGlobalTime = track->GetGlobalTime() / CLHEP::microsecond;
    const G4ThreeVector& momentum = track->GetMomentumDirection();

    AddHit(hitPosition, energy, hitGlobalTime);  // this increases fNHits

    fProcessID.emplace_back(processID);
    fVolumeID.emplace_back(geometryInfo.GetIDFromVolume(volumeName));
    fKineticEnergy.emplace_back(track->GetKineticEnergy() / CLHEP::keV);
    fMomentumDirection.emplace_back(momentum.x(), momentum.y(), momentum.z());

    string isotopeName;
    int atomicNumber = 0;
    int atomicMassNumber = 0;

    if (metadata->GetStoreHadronicTargetInfo() && track->GetCurrentStepNumber() != 0 &&
        process->GetProcessType() == G4ProcessType::fHadronic) {
        auto hadronicProcess = dynamic_cast<const G4HadronicProcess*>(process);
        G4Nucleus nucleus = *(hadronicProcess->GetTargetNucleus());

        auto isotope = nucleus.GetIsotope();
        if (isotope) {
            isotopeName = isotope->GetName();
            atomicNumber = isotope->GetZ();
            atomicMassNumber = isotope->GetN();
        }
    }

    if (metadata->GetStoreHadronicTargetInfo()) {
        fHadronicTargetIsotopeName.emplace_back(isotopeName);
        fHadronicTargetIsotopeZ.emplace_back(atomicNumber);
        fHadronicTargetIsotopeA.emplace_back(atomicMassNumber);
    }

    SimulationManager::GetOutputManager()->AddEnergyToVolumeForParticleForProcess(energy, volumeName,
                                                                                  particleName, processName);
}

void OutputManager::RemoveUnwantedTracks() {
    const auto& metadata = fSimulationManager->GetRestMetadata();
    set<int> trackIDsToKeep;  // We populate this container with the tracks we want to keep
    for (const auto& track : fEvent->fTracks) {
        // If one children track is kept, we keep all the parents
        if (trackIDsToKeep.count(track.GetTrackID()) > 0) {
            continue;
        }
        const auto hits = track.GetHits();
        for (int i = 0; i < int(hits.GetNumberOfHits()); i++) {
            const auto energy = hits.GetEnergy(i);
            if (!fSimulationManager->GetRestMetadata()->GetRemoveUnwantedTracksKeepZeroEnergyTracks() &&
                energy <= 0) {
                continue;
            }
            const auto volume = metadata->GetGeant4GeometryInfo().GetVolumeFromID(hits.GetVolumeId(i));
            if (metadata->IsKeepTracksVolume(volume)) {
                trackIDsToKeep.insert(track.GetTrackID());
                auto parentTrack = track.GetParentTrack();
                while (parentTrack != nullptr) {
                    trackIDsToKeep.insert(parentTrack->GetTrackID());
                    parentTrack = parentTrack->GetParentTrack();
                }
            }
        }
    }
    // const size_t numberOfTracksBefore = fEvent->fTracks.size();

    vector<TRestGeant4Track> tracksAfterRemoval;
    for (const auto& track : fEvent->fTracks) {
        // we do this to preserve original order
        if (trackIDsToKeep.count(track.GetTrackID()) > 0) {
            tracksAfterRemoval.push_back(*(fEvent->GetTrackByID(track.GetTrackID())));
        }
    }

    fEvent->fTracks = tracksAfterRemoval;

    // Updated indices
    fEvent->fTrackIDToTrackIndex.clear();
    for (int i = 0; i < int(fEvent->fTracks.size()); i++) {
        fEvent->fTrackIDToTrackIndex[fEvent->fTracks[i].GetTrackID()] = i;
    }

    /*
    const size_t numberOfTracksAfter = fEvent->fTracks.size();
    cout << "EventID: " << fEvent->GetID() << " Removed " << numberOfTracksBefore - numberOfTracksAfter
         << " tracks out of " << numberOfTracksBefore << endl;
     */
}
