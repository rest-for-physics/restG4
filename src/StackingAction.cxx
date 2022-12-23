
#include "StackingAction.h"

#include <G4ParticleTable.hh>
#include <G4ParticleTypes.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>
#include <G4VProcess.hh>

#include "SimulationManager.h"

StackingAction::StackingAction(SimulationManager* simulationManager) : fSimulationManager(simulationManager) {
    fMaxAllowedLifetime = 100;  // TODO: update this

    fMaxAllowedLifetimeWithUnit = G4BestUnit(fMaxAllowedLifetime, "Time");

    fParticlesToIgnore = {
        G4NeutrinoE::Definition(),      G4AntiNeutrinoE::Definition(), G4NeutrinoMu::Definition(),
        G4AntiNeutrinoMu::Definition(), G4NeutrinoTau::Definition(),   G4AntiNeutrinoTau::Definition(),
    };

    /*
for (const auto& particle : fParticlesToIgnore) {
}
    */
}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track) {
    const G4ClassificationOfNewTrack decayClassification =
        fSimulationManager->GetRestMetadata()->isFullChainActivated() ? fWaiting : fKill;
    if (track->GetParentID() <= 0) {
        // always process the first track regardless
        return fUrgent;
    }

    auto particle = track->GetParticleDefinition();
    if (fParticlesToIgnore.find(particle) != fParticlesToIgnore.end()) {
        // ignore this track
        return fKill;
    }

    if (track->GetCreatorProcess()->GetProcessType() != G4ProcessType::fDecay) {
        return fUrgent;
    }

    if (particle->GetParticleType() == "nucleus" && !particle->GetPDGStable()) {
        // unstable nucleus
        if (particle->GetPDGLifeTime() > fMaxAllowedLifetime) {
            G4String energy = G4BestUnit(track->GetKineticEnergy(), "Energy");
            G4String lifeTime = G4BestUnit(particle->GetPDGLifeTime(), "Time");
            return decayClassification;
        }
    }

    return fUrgent;
}

void StackingAction::NewStage() {
    /*
     * Called after processing all Urgent tracks
     * Close event and start a new sub event if there are waiting tracks
     */

    const auto outputManager = fSimulationManager->GetOutputManager();

    const Int_t subEventID = outputManager->fEvent->GetSubID();
    outputManager->FinishAndSubmitEvent();
    outputManager->fEvent->SetSubID(subEventID + 1);
}

StackingAction::~StackingAction() = default;
