
#ifndef REST_SIMULATIONMANAGER_H
#define REST_SIMULATIONMANAGER_H

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

class SimulationManager {
   public:
    SimulationManager();
    ~SimulationManager();

    TRestRun* fRestRun = nullptr;
    TRestGeant4Track* fRestGeant4Track = nullptr;
    TRestGeant4Event *fRestGeant4Event = nullptr, *fRestGeant4SubEvent = nullptr;
    TRestGeant4Metadata* fRestGeant4Metadata = nullptr;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists = nullptr;

    TH1D initialEnergySpectrum;
    TH1D initialAngularDistribution;
};

#endif  // REST_SIMULATIONMANAGER_H
