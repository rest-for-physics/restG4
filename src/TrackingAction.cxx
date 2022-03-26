

#include "TrackingAction.h"

#include <G4ParticleTypes.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>

#include "EventAction.h"
#include "RunAction.h"

extern TRestGeant4Metadata* restG4Metadata;
extern TRestGeant4Event* restG4Event;
extern TRestGeant4Track* restTrack;

G4double prevTime = 0;
G4String aux;

TrackingAction::TrackingAction(RunAction* RA, EventAction* EA)
    : G4UserTrackingAction(),
      fRun(RA),
      fEvent(EA)

{
    fFullChain = false;

    fFullChain = restG4Metadata->isFullChainActivated();

    if (fFullChain)
        G4cout << "Full chain is active" << G4endl;
    else
        G4cout << "Full chain is NOT active" << G4endl;
}

TrackingAction::~TrackingAction() {}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    if (restG4Metadata->GetVerboseLevel() >= REST_Extreme)
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

    Int_t ID = track->GetTrackID();
    if (!fFullChain && fCharge > 2 && ID > 1 && G4StrUtil::contains(name, "[")) {
        auto tr = (G4Track*)track;
        tr->SetTrackStatus(fStopAndKill);
    }
}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
    restTrack->SetTrackTimeLength(track->GetLocalTime() / microsecond);

    //   G4cout << "Storing track : Number of hits : " <<
    //   restTrack->GetNumberOfHits() << G4endl;

    // if (restTrack->GetNumberOfHits() > 0 || restG4Metadata->RegisterEmptyTracks())
    restG4Event->AddTrack(*restTrack);
}
