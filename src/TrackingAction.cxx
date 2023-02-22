

#include "TrackingAction.h"

#include <G4RunManager.hh>
#include <G4Track.hh>

#include "SimulationManager.h"

using namespace std;

TrackingAction::TrackingAction(SimulationManager *simulationManager)
        : G4UserTrackingAction(), fSimulationManager(simulationManager) {}

void TrackingAction::PreUserTrackingAction(const G4Track *track) {
    // G4cout << track->GetVolume()->GetLogicalVolume()->GetName() << G4endl;
    fSimulationManager->GetOutputManager()->RecordTrack(track);
}

void TrackingAction::PostUserTrackingAction(const G4Track *track) {
    fSimulationManager->GetOutputManager()->UpdateTrack(track);
}
