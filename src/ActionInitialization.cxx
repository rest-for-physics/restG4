//
// Created by lobis on 6/14/2022.
//

#include "ActionInitialization.h"

#include "DetectorConstruction.h"
#include "EventAction.h"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "SimulationManager.h"
#include "StackingAction.h"
#include "SteppingAction.h"
#include "SteppingVerbose.h"
#include "TrackingAction.h"

using namespace std;

ActionInitialization::ActionInitialization(SimulationManager* simulationManager)
    : G4VUserActionInitialization(), fSimulationManager(simulationManager) {}

ActionInitialization::~ActionInitialization() = default;

void ActionInitialization::BuildForMaster() const { SetUserAction(new RunAction(fSimulationManager)); }

void ActionInitialization::Build() const {
    fSimulationManager->InitializeOutputManager();

    SetUserAction(new PrimaryGeneratorAction(fSimulationManager));
    SetUserAction(new RunAction(fSimulationManager));
    SetUserAction(new EventAction(fSimulationManager));
    SetUserAction(new SteppingAction(fSimulationManager));
    SetUserAction(new StackingAction(fSimulationManager));
    SetUserAction(new TrackingAction(fSimulationManager));

    /*
    G4EventManager::GetEventManager()->SetNumberOfAdditionalWaitingStacks(1);  // optical stack
     */
}

G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const {
    return new SteppingVerbose(fSimulationManager);
}