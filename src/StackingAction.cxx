
#include "StackingAction.h"

#include <G4ParticleTable.hh>
#include <G4ParticleTypes.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>
#include <G4VProcess.hh>
#include <regex>

#include "SimulationManager.h"

using namespace std;

StackingAction::StackingAction(SimulationManager* simulationManager) : fSimulationManager(simulationManager) {
    fMaxAllowedLifetime = fSimulationManager->GetRestMetadata()->GetSubEventTimeDelay() * CLHEP::microsecond;

    fMaxAllowedLifetimeWithUnit = G4BestUnit(fMaxAllowedLifetime, "Time");

    fParticlesToIgnore = {
        G4NeutrinoE::Definition(),      G4AntiNeutrinoE::Definition(), G4NeutrinoMu::Definition(),
        G4AntiNeutrinoMu::Definition(), G4NeutrinoTau::Definition(),   G4AntiNeutrinoTau::Definition(),
    };
}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track) {
    const bool isFullChain = fSimulationManager->GetRestMetadata()->isFullChainActivated();
    const G4ClassificationOfNewTrack decayClassification = isFullChain ? fWaiting : fKill;
    if (track->GetParentID() <= 0) {
        // always process the first track regardless
        return fUrgent;
    }

    auto particle = track->GetParticleDefinition();
    if (fParticlesToIgnore.find(particle) != fParticlesToIgnore.end()) {
        // ignore this track
        return fKill;
    }

    if (isFullChain) {
        const string particleName = particle->GetParticleName();
        // the particle may be an excited state of an isotope we want to remove
        // strip the excited state which is assumed to always be inside some "[]"
        regex pattern("\\[.*\\]");
        std::string particleNameStripped = regex_replace(particleName, pattern, "");
        if (fSimulationManager->GetRestMetadata()->IsIsotopeFullChainStop(particleNameStripped)) {
            return fKill;
        }
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

    const auto outputManager = SimulationManager::GetOutputManager();

    const Int_t subEventID = outputManager->fEvent->GetSubID();
    outputManager->FinishAndSubmitEvent();
    outputManager->fEvent->SetSubID(subEventID + 1);
}

StackingAction::~StackingAction() = default;
