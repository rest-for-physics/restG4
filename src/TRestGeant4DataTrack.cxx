//
// Created by lobis on 12/10/2021.
//

#include <SteppingAction.h>
#include <TRestGeant4DataEvent.h>
#include <TRestGeant4DataSteps.h>
#include <TRestGeant4DataTrack.h>
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

using namespace std;

TRestGeant4DataTrack::TRestGeant4DataTrack(const G4Track* track) {
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
    }

    fInitialKineticEnergy = track->GetKineticEnergy() / CLHEP::keV;

    fWeight = track->GetWeight();

    //
    G4String energyWithUnits = G4BestUnit(fInitialKineticEnergy * CLHEP::keV, "Energy");

    spdlog::debug(
        "DataModelTrack::DataModelTrack - Track ID {} - Parent ID {} - Particle {} ({}) created by {} - "
        "Energy {}",
        fTrackID, fParentID, fParticleName, fParticleID,
        (fCreatorProcess.IsNull() ? "IS-PRIMARY-PARTICLE" : fCreatorProcess), energyWithUnits);
}

void TRestGeant4DataTrack::InsertStep(const G4Step* step) { fSteps.InsertStep(step); }

void TRestGeant4DataTrack::UpdateTrack(const G4Track* track) {
    if (track->GetTrackID() != fTrackID) {
        spdlog::error("DataModelTrack::UpdateTrack - mismatch of trackID!");
        exit(1);
    }

    fTrackLength = track->GetTrackLength() / CLHEP::mm;

    auto steppingAction = (SteppingAction*)G4EventManager::GetEventManager()->GetUserSteppingAction();
    auto secondaries = steppingAction->GetfSecondary();
    fNumberOfSecondaries = (int)secondaries->size();
}
