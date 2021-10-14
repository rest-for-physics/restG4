
#include "SteppingAction.h"

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4Track.h>

#include <G4DynamicParticle.hh>
#include <G4OpticalPhoton.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4SteppingManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

#include "GlobalManager.h"
#include "OutputManager.h"

extern Int_t biasing;

SteppingAction::SteppingAction()
    : fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()),
      fOutputManager(OutputManager::Instance()) {
    if (biasing > 1) restBiasingVolume = fRestGeant4Metadata->GetBiasingVolume(biasing - 1);
}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto track = step->GetTrack();

    // optical photons
    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        // TODO: Do something with optical photons
    }

    fOutputManager->RecordStep(step);

    return;
}

G4TrackVector* SteppingAction::GetfSecondary() {
    G4SteppingManager* steppingManager = fpSteppingManager;
    return steppingManager->GetfSecondary();
}
