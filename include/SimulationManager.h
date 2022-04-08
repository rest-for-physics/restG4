
#ifndef REST_SIMULATIONMANAGER_H
#define REST_SIMULATIONMANAGER_H

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

class SimulationManager {
   public:
    static SimulationManager* Instance();
    ~SimulationManager();

    TRestRun* fRestRun;
    TRestGeant4Track* fRestGeant4Track;
    TRestGeant4Event *fRestGeant4Event, *fRestGeant4SubEvent;
    TRestGeant4Metadata* fRestGeant4Metadata;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists;
    Int_t fBiasing = 0;

   private:
    SimulationManager();
    static SimulationManager* fInstance;
};

#endif  // REST_SIMULATIONMANAGER_H
