//
// Created by lobis on 12/10/2021.
//

#include "SensitiveDetector.h"

#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>
#include <string>

#include "OutputManager.h"
#include "spdlog/spdlog.h"

using namespace std;

SensitiveDetector::SensitiveDetector(const G4String& name)
    : G4VSensitiveDetector(name), fOutputManager(OutputManager::Instance()) {}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    auto energy = step->GetTotalEnergyDeposit() / keV;
    fOutputManager->AddSensitiveEnergy(energy);
    return true;  // return value will always be ignored (legacy)
}