
#ifndef REST_SIMULATIONMANAGER_H
#define REST_SIMULATIONMANAGER_H

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

#include <queue>

class OutputManager;

class SimulationManager {
   public:
    SimulationManager();
    ~SimulationManager();

    TRestRun* fRestRun = nullptr;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists = nullptr;
    TRestGeant4Metadata* fRestGeant4Metadata = nullptr;

    TH1D initialEnergySpectrum;
    TH1D initialAngularDistribution;

    void InitializeOutputManager();
    static OutputManager* GetOutputManager() { return fOutputManager; }

    TRestGeant4Event fEvent;  // Branch on EventTree

   private:
    static thread_local OutputManager* fOutputManager;
    std::queue<std::unique_ptr<TRestGeant4Event> > fEventContainer;
};

class OutputManager {
   public:
    OutputManager(const SimulationManager*);
    void UpdateEvent();

   private:
    std::unique_ptr<TRestGeant4Event> fEvent{};
    const SimulationManager* fSimulationManager = nullptr;
};

#endif  // REST_SIMULATIONMANAGER_H
