//
// Created by lobis on 12/10/2021.
//

#include "StackingAction.h"

#include <GlobalManager.h>

#include <G4ParticleTable.hh>
#include <G4ParticleTypes.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>
#include <G4VProcess.hh>

#include "OutputManager.h"
#include "spdlog/spdlog.h"

StackingAction::StackingAction()
    : G4UserStackingAction(),
      fOutputManager(OutputManager::Instance()),
      fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()),
      fMaxAllowedLifetime(100 / nanosecond) {
    fParticlesToIgnore = {
        G4NeutrinoE::Definition(),      G4AntiNeutrinoE::Definition(), G4NeutrinoMu::Definition(),
        G4AntiNeutrinoMu::Definition(), G4NeutrinoTau::Definition(),   G4AntiNeutrinoTau::Definition(),
    };

    fMaxAllowedLifetime = fRestGeant4Metadata->GetSubEventTimeDelay() / microsecond;
    fFullChain = fRestGeant4Metadata->isFullChainActivated();
    // Debugging
    if (!G4Threading::IsMasterThread() || !G4Threading::IsMultithreadedApplication()) {
        for (const auto& particle : fParticlesToIgnore) {
            spdlog::info("StackingAction::StackingAction ---> Particle to ignore: {}",
                         particle->GetParticleName());
        }

        spdlog::info("StackingAction::StackingAction - 'fFullChain' is set to {}", fFullChain);

        fMaxAllowedLifetimeWithUnit = G4BestUnit(fMaxAllowedLifetime, "Time");
        spdlog::info("StackingAction::StackingAction - 'fMaxAllowedLifetime' is {}",
                     fMaxAllowedLifetimeWithUnit);
    }
}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track) {
    if (track->GetParentID() <= 0) {
        // always process the first track regardless
        return fUrgent;
    }

    auto particle = track->GetParticleDefinition();
    if (fParticlesToIgnore.find(particle) != fParticlesToIgnore.end()) {
        // ignore this track
        spdlog::debug(
            "StackingAction::ClassifyNewTrack ---> Killing track with '{}' (is in 'fParticlesToIgnore')",
            particle->GetParticleName());
        return fKill;
    }

    if (particle == G4OpticalPhoton::Definition()) {
        // return fKill;
    }

    if (track->GetCreatorProcess()->GetProcessType() != G4ProcessType::fDecay) {
        return fUrgent;
    }

    if (particle->GetParticleType() == "nucleus" && !particle->GetPDGStable()) {
        // unstable nucleus
        if (particle->GetPDGLifeTime() > fMaxAllowedLifetime) {
            G4String energy = G4BestUnit(track->GetKineticEnergy(), "Energy");
            G4String lifeTime = G4BestUnit(particle->GetPDGLifeTime(), "Time");

            spdlog::debug(
                "StackingAction::ClassifyNewTrack ---> Splitting unstable '{}' track (KE: {}, lifetime of {} "
                "> {})",
                particle->GetParticleName(), energy, lifeTime, fMaxAllowedLifetimeWithUnit);

            if (track->GetKineticEnergy() > 0) {
                spdlog::debug(
                    "StackingAction::ClassifyNewTrack ---> Killing unstable '{}' track (KE: {}). It's energy "
                    "is not zero.",
                    particle->GetParticleName(), energy);
            }

            //  TODO: kill/wait it when it stops, not before
            return fFullChain ? fWaiting : fKill;
        }
    }
    return fUrgent;
}

void StackingAction::NewStage() {
    /*
     * Called after processing all Urgent tracks
     * Close event and start a new sub event if there are waiting tracks
     */
    spdlog::debug("StackingAction::NewStage");

    const Int_t subEventID = fOutputManager->GetSubEventID();
    fOutputManager->FinishAndSubmitEvent();
    fOutputManager->SetSubEventID(subEventID + 1);
}

StackingAction::~StackingAction() = default;