//
// Created by lobis on 12/10/2021.
//

#include <SteppingAction.h>
#include <TRestGeant4Event.h>
#include <TRestGeant4Hits.h>
#include <TRestGeant4Track.h>
#include <spdlog/spdlog.h>

#include <G4Event.hh>
#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include <G4ParticleDefinition.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4Scintillation.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4Transportation.hh>
#include <G4UnitsTable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VProcess.hh>
#include <G4VUserPrimaryVertexInformation.hh>

#include "GlobalManager.h"
#include "OutputManager.h"

using namespace std;

TRestGeant4Track::TRestGeant4Track(const G4Track* track) {
    fTrackID = track->GetTrackID();
    fParentID = track->GetParentID();

    auto particle = track->GetParticleDefinition();
    fParticleName = particle->GetParticleName();
    fParticleID = particle->GetPDGEncoding();
    fParticleType = particle->GetParticleType();
    fParticleSubType = particle->GetParticleSubType();

    auto creatorProcess = track->GetCreatorProcess();
    if (creatorProcess) {
        fCreatorProcess = creatorProcess->GetProcessName();
    } else {
        fCreatorProcess = "IS-PRIMARY-PARTICLE";
    }

    fInitialKineticEnergy = track->GetKineticEnergy() / CLHEP::keV;

    fWeight = track->GetWeight();

    fSecondaryTrackIDs = {};  // To be filled at TRestGeant4Track::UpdateTrack

    G4String energyWithUnits = G4BestUnit(fInitialKineticEnergy * CLHEP::keV, "Energy");

    spdlog::debug(
        "DataModelTrack::DataModelTrack - Track ID {} - Parent ID {} - Particle {} ({}) created by {} - "
        "Energy {}",
        fTrackID, fParentID, fParticleName, fParticleID,
        (fCreatorProcess.IsNull() ? "IS-PRIMARY-PARTICLE" : fCreatorProcess), energyWithUnits);
}

void TRestGeant4Track::InsertStep(const G4Step* step) { fHits.InsertStep(step); }

void TRestGeant4Track::UpdateTrack(const G4Track* track) {
    if (track->GetTrackID() != fTrackID) {
        spdlog::error("DataModelTrack::UpdateTrack - mismatch of trackID!");
        exit(1);
    }

    fTrackLength = track->GetTrackLength() / CLHEP::mm;

    // auto steppingAction = (SteppingAction*)G4EventManager::GetEventManager()->GetUserSteppingAction();
    // auto secondaries = steppingAction->GetSecondary();
}

void TRestGeant4Track::AddSecondary(Int_t trackID) {
    // This method is called only for recorded events, at the end of event action
    for (const auto& existingTrackID : fSecondaryTrackIDs) {
        if (existingTrackID == trackID) {
            spdlog::error("TRestGeant4Track::AddSecondary - Adding the same secondary twice, please check!");
            exit(1);
        }
    }

    fSecondaryTrackIDs.push_back(trackID);
}