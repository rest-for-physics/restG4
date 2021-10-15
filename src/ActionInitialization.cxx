//
// Created by lobis on 12/10/2021.
//

#include "ActionInitialization.h"

#include "DetectorConstruction.h"
#include "EventAction.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "StackingAction.h"
#include "SteppingAction.h"
#include "SteppingVerbose.h"
#include "TrackingAction.h"
#include "spdlog/spdlog.h"

ActionInitialization::ActionInitialization() : G4VUserActionInitialization() {
    spdlog::debug("ActionInitialization::ActionInitialization");
}

ActionInitialization::~ActionInitialization() = default;

void ActionInitialization::BuildForMaster() const {
    spdlog::info("ActionInitialization::BuildForMaster");
    SetUserAction(new RunAction);
}

void ActionInitialization::Build() const {
    spdlog::info("ActionInitialization::Build");

    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
    SetUserAction(new SteppingAction);
    // SetUserAction(new StackingAction);
    SetUserAction(new TrackingAction);

    // G4EventManager::GetEventManager()->SetNumberOfAdditionalWaitingStacks(1);  // optical stack
}

G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const { return new SteppingVerbose(); }
