

#include "TrackingAction.h"

#include <G4ParticleTypes.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>

#include "SimulationManager.h"

using namespace std;

TrackingAction::TrackingAction(SimulationManager* simulationManager)
    : G4UserTrackingAction(), fSimulationManager(simulationManager) {}

TrackingAction::~TrackingAction() {}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    fSimulationManager->GetOutputManager()->RecordTrack(track);
}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
    fSimulationManager->GetOutputManager()->UpdateTrack(track);
}
