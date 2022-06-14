

#include "TrackingAction.h"

#include <G4ParticleTypes.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>

#include "EventAction.h"
#include "RunAction.h"
#include "SimulationManager.h"

using namespace std;

G4double prevTime = 0;
G4String aux;

TrackingAction::TrackingAction(SimulationManager* simulationManager, RunAction* runAction,
                               EventAction* eventAction)
    : G4UserTrackingAction(),
      fSimulationManager(simulationManager),
      fRun(runAction),
      fEvent(eventAction)

{
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;

    fFullChain = false;

    fFullChain = restG4Metadata->isFullChainActivated();

    if (fFullChain)
        G4cout << "Full chain is active" << G4endl;
    else
        G4cout << "Full chain is NOT active" << G4endl;
}

TrackingAction::~TrackingAction() {}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    TRestRun* restRun = fSimulationManager->fRestRun;
    TRestGeant4Track* restTrack = fSimulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = fSimulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = fSimulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;
    TRestGeant4PhysicsLists* restPhysList = fSimulationManager->fRestGeant4PhysicsLists;
    Int_t& biasing = fSimulationManager->fBiasing;

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) {
        if (track->GetTrackID() % 10 == 0) {
            cout << "EXTREME: Processing track " << track->GetTrackID() << endl;
        }
    }
    G4ParticleDefinition* particle = track->GetDefinition();
    const auto& name = particle->GetParticleName();
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
    if (!fFullChain && fCharge > 2 && ID > 1 &&
#ifdef GEANT4_VERSION_LESS_11_0_0
        !name.contains("[")
#else
        !G4StrUtil::contains(name, "[")
#endif
    ) {
        auto tr = (G4Track*)track;
        tr->SetTrackStatus(fStopAndKill);
    }
}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
    TRestGeant4Track* restTrack = fSimulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = fSimulationManager->fRestGeant4Event;

    restTrack->SetTrackTimeLength(track->GetLocalTime() / microsecond);

    //   G4cout << "Storing track : Number of hits : " <<
    //   restTrack->GetNumberOfHits() << G4endl;

    // if (restTrack->GetNumberOfHits() > 0 || restG4Metadata->RegisterEmptyTracks())
    restG4Event->AddTrack(*restTrack);
}
