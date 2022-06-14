
#include "SimulationManager.h"

#include <G4Threading.hh>

using namespace std;

SimulationManager::SimulationManager() {
    // Only master thread should create the SimulationManager
    if (!G4Threading::IsMasterThread()) {
        cout << "Only master thread should create the SimulationManager!" << endl;
        exit(1);
    }
}

SimulationManager::~SimulationManager() = default;
