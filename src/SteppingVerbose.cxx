//
// Created by lobis on 6/20/2022.
//

#include "SteppingVerbose.h"

#include <SimulationManager.h>

#include <G4OpticalPhoton.hh>
#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4SteppingManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>

#include "DetectorConstruction.h"
#include "EventAction.h"

SteppingVerbose::SteppingVerbose(SimulationManager* simulationManager)
    : G4SteppingVerbose(), fSimulationManager(simulationManager) {}
SteppingVerbose::~SteppingVerbose() {}

void SteppingVerbose::TrackingStarted() {
    CopyState();
    fSimulationManager->GetOutputManager()->RecordStep(fStep);
}

void SteppingVerbose::StepInfo() {}
