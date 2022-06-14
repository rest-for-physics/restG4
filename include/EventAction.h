
#ifndef EventAction_h
#define EventAction_h 1

#include <TH1D.h>
#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4Track.h>
#include <TStopwatch.h>

#include <G4UserEventAction.hh>
#include <globals.hh>

class SimulationManager;

class EventAction : public G4UserEventAction {
   public:
    EventAction(SimulationManager*);
    ~EventAction();

   public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

   private:
    SimulationManager* fSimulationManager;
    TStopwatch fTimer;

    Double_t absDouble(Double_t x) {
        if (x > 0) return x;
        return -x;
    }

    int SetTrackSubEventIDs();
    static void FillSubEvent(Int_t subId);
    // old method `FillSubEvent` has been split into `FillSubEvent` and `ReOrderTrackIds` for speed
    static void ReOrderTrackIds(Int_t subId);

    // variable used to track the number of events that hit the sensitive volume
    UInt_t sensitive_volume_hits_count = 0;
};

#endif
