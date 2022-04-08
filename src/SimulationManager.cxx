
#include "SimulationManager.h"

#include <G4Threading.hh>

using namespace std;

SimulationManager* SimulationManager::fInstance = nullptr;

SimulationManager* SimulationManager::Instance() {
    if (!fInstance) {
        fInstance = new SimulationManager();
    }
    return fInstance;
}

SimulationManager::SimulationManager() {
    // Master thread should create the SimulationManager, worker threads should spawn after the manager is
    // created
    if (!G4Threading::IsMasterThread()) {
        exit(1);
    }
}

SimulationManager::~SimulationManager() = default;
