
#ifndef TrackingAction_h
#define TrackingAction_h 1

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4Track.h>

#include <G4Track.hh>
#include <G4UserTrackingAction.hh>
#include <globals.hh>

class RunAction;
class EventAction;
class SimulationManager;

class TrackingAction : public G4UserTrackingAction {
   public:
    TrackingAction(SimulationManager*);

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

   private:
    SimulationManager* fSimulationManager;
};

#endif
