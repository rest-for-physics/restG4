
#include "TrackingAction.h"

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

extern TRestGeant4Event* restG4Event;
extern TRestGeant4Track* restTrack;

TrackingAction::TrackingAction()
    : G4UserTrackingAction(),
      fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()),
      fOutputManager(OutputManager::Instance())

{
    fRun = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
    fEvent = (EventAction*)G4RunManager::GetRunManager()->GetUserEventAction();

    fFullChain = false;

    fFullChain = fRestGeant4Metadata->isFullChainActivated();

    if (fFullChain)
        G4cout << "Full chain is active" << G4endl;
    else
        G4cout << "Full chain is NOT active" << G4endl;
}

TrackingAction::~TrackingAction() = default;

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    fOutputManager->RecordTrack(track);
    return;

    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Extreme)
        if (track->GetTrackID() % 10 == 0) {
            cout << "EXTREME: Processing track " << track->GetTrackID() << endl;
        }
    G4ParticleDefinition* particle = track->GetDefinition();
    G4String name = particle->GetParticleName();
    fCharge = particle->GetPDGCharge();

    restTrack->RemoveHits();

    restTrack->SetTrackID(track->GetTrackID());
    restTrack->SetParentID(track->GetParentID());
    restTrack->SetKineticEnergy(track->GetKineticEnergy() / keV);
    restTrack->SetParticleName(name);
    restTrack->SetGlobalTrackTime(track->GetGlobalTime() / second);

    G4ThreeVector trkOrigin = track->GetPosition();
    restTrack->SetTrackOrigin(trkOrigin.x(), trkOrigin.y(), trkOrigin.z());

    // We finish after the de-excitation of the resulting nucleus (we skip the
    // full chain, just first decay) On future we must add an option through
    // TRestGeant4Metadata to store a given number of decays

    //   if( fFullChain == true ) G4cout << "Full chain active" << G4endl;
    //   else  G4cout << "Full chain not active" << G4endl;

    Int_t ID = track->GetTrackID();
    if (!fFullChain && fCharge > 2 && ID > 1 && !name.contains("[")) {
        auto tr = (G4Track*)track;
        tr->SetTrackStatus(fStopAndKill);
    }

    /*
    if ( fFullChain == true && fCharge > 2  && ID > 1 && !name.contains("["))
    {
        restTrack->IncreaseSubEventID();
    }
    */
}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
    auto particle = track->GetParticleDefinition();

    auto steppingAction = (SteppingAction*)G4EventManager::GetEventManager()->GetUserSteppingAction();
    auto fSecondary = steppingAction->GetfSecondary();

    /*
    spdlog::debug("TrackingAction::PostUserTrackingAction ---> TrackID: {:4d}, ParentID: {:4d}, Particle:
    {:>10}, Secondaries: {:2d}", track->GetTrackID(), track->GetParentID(), particle->GetParticleName(),
    fSecondary->size());
    */

    for (const auto& secondary : *fSecondary) {
        G4String energyWithUnits = G4BestUnit(secondary->GetKineticEnergy(), "Energy");
        auto creatorProcess = secondary->GetCreatorProcess();
        /*
        spdlog::debug("TrackingAction::PostUserTrackingAction ---> ---> Secondary: {:>10}, KE:  {:>15},
        Creator Process: {:>20} ({})", secondary->GetDynamicParticle()->GetDefinition()->GetParticleName(),
        energyWithUnits, creatorProcess->GetProcessName(),
                      G4VProcess::GetProcessTypeName(creatorProcess->GetProcessType()));
                      */
    }

    fOutputManager->UpdateTrack(track);

    return;
    restTrack->SetTrackTimeLength(track->GetLocalTime() / microsecond);

    //   G4cout << "Storing track : Number of hits : " <<
    //   restTrack->GetNumberOfHits() << G4endl;

    // if (restTrack->GetNumberOfHits() > 0 || fRestGeant4Metadata->RegisterEmptyTracks())
    restG4Event->AddTrack(*restTrack);
}
