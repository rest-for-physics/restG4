
#include "SensitiveDetector.h"

#include <G4Geantino.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>

#include "SimulationManager.h"

using namespace std;

SensitiveDetector::SensitiveDetector(SimulationManager* simulationManager, const G4String& name)
    : fSimulationManager(simulationManager), G4VSensitiveDetector(name) {}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    // return value will always be ignored, its present for backwards compatibility (I guess)
    const bool isGeantino = step->GetTrack()->GetParticleDefinition() == G4Geantino::Definition();

    if (isGeantino) {
        const auto volumeName = (TString)step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
        const auto length = step->GetStepLength() / CLHEP::mm;
        fSimulationManager->GetOutputManager()->AddSensitiveEnergy(length, volumeName);
        return true;
    } else {
        auto energy = step->GetTotalEnergyDeposit() / keV;

        if (energy <= 0) {
            return true;
        }

        const auto volumeName = (TString)step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

        fSimulationManager->GetOutputManager()->AddSensitiveEnergy(energy, volumeName);

        return true;
    }
}