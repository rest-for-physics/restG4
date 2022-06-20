
#include "SimulationManager.h"

#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4Threading.hh>

using namespace std;

SimulationManager::SimulationManager() {
    // Only master thread should create the SimulationManager
    if (!G4Threading::IsMasterThread()) {
        cout << "Only master thread should create the SimulationManager!" << endl;
        exit(1);
    }
}

SimulationManager::~SimulationManager() {
    delete fRestRun;
    delete fRestGeant4Metadata;
    delete fRestGeant4PhysicsLists;
}

// OutputManager
void OutputManager::UpdateEvent() {
    auto event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    fEvent = make_unique<TRestGeant4Event>(event, *fSimulationManager->fRestGeant4Metadata);
}

// Geant4Lib constructors
TRestGeant4Event::TRestGeant4Event(const G4Event* event, const TRestGeant4Metadata& metadata)
    : TRestGeant4Event() {
    SetID(event->GetEventID());
    SetOK(true);
    time_t system_time = time(nullptr);

    SetTime((Double_t)system_time);

    // Defining if the hits in a given volume will be stored
    for (int i = 0; i < metadata.GetNumberOfActiveVolumes(); i++) {
        if (metadata.GetStorageChance(i) >= 1.00) {
            ActivateVolumeForStorage(i);
        } else {
            Double_t randomNumber = G4UniformRand();
            if (metadata.GetStorageChance(i) >= randomNumber) {
                ActivateVolumeForStorage(i);
            } else {
                DisableVolumeForStorage(i);
            }
        }
    }
}