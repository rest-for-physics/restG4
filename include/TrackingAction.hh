
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

class TrackingAction : public G4UserTrackingAction {
   public:
    TrackingAction(RunAction*, EventAction*);
    ~TrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

   private:
    RunAction* fRun;
    EventAction* fEvent;

    G4double fCharge, fMass;

    G4bool fFullChain;

    Double_t fGlobalTime;
};

#endif
