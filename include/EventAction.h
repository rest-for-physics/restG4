
#ifndef EventAction_h
#define EventAction_h 1

using namespace std;

#include <TH1D.h>
#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4Track.h>

#include <G4UserEventAction.hh>
#include <globals.hh>

class OutputManager;

class EventAction : public G4UserEventAction {
   public:
    EventAction();
    ~EventAction();

   public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

   private:
    TRestGeant4Metadata* fRestGeant4Metadata;
    OutputManager* fOutputManager;

    int SetTrackSubEventIDs();
    void FillSubEvent(Int_t subId);
    // old method `FillSubEvent` has been split into `FillSubEvent` and `ReOrderTrackIds` for speed
    void ReOrderTrackIds(Int_t subId);

    // variable used to track the number of events that hit the sensitive volume
    UInt_t sensitive_volume_hits_count = 0;
};

#endif
