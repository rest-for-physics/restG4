
#include "TrackingAction.h"

#include <spdlog/spdlog.h>

#include <G4ParticleTypes.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>

#include "EventAction.h"
#include "GlobalManager.h"
#include "OutputManager.h"
#include "RunAction.h"
#include "SteppingAction.h"

TrackingAction::TrackingAction()
    : G4UserTrackingAction(),
      fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()),
      fOutputManager(OutputManager::Instance())

{
    return;
}

TrackingAction::~TrackingAction() = default;

void TrackingAction::PreUserTrackingAction(const G4Track* track) { fOutputManager->RecordTrack(track); }

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
    auto particle = track->GetParticleDefinition();

    auto steppingAction = (SteppingAction*)G4EventManager::GetEventManager()->GetUserSteppingAction();
    auto secondaries = steppingAction->GetSecondary();

    spdlog::debug(
        "TrackingAction::PostUserTrackingAction ---> TrackID: {:4d}, ParentID: {:4d}, Particle: {:>10}, "
        "Secondaries: {:2d}",
        track->GetTrackID(), track->GetParentID(), particle->GetParticleName(), secondaries->size());

    for (const auto& secondary : *secondaries) {
        G4String energyWithUnits = G4BestUnit(secondary->GetKineticEnergy(), "Energy");
        auto creatorProcess = secondary->GetCreatorProcess();

        spdlog::debug(
            "TrackingAction::PostUserTrackingAction ---> ---> Secondary ID: {:>5} Particle:{:>10}, KE: "
            "{:>15}, Creator Process: {:>20} ({})",
            secondary->GetTrackID(), secondary->GetDynamicParticle()->GetDefinition()->GetParticleName(),
            energyWithUnits, creatorProcess->GetProcessName(),
            G4VProcess::GetProcessTypeName(creatorProcess->GetProcessType()));
    }

    fOutputManager->UpdateTrack(track);
}
