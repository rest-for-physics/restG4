
#include "SensitiveDetector.h"

#include <G4Geantino.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>

#include "SimulationManager.h"

using namespace std;

SensitiveDetector::SensitiveDetector(SimulationManager* simulationManager, const G4String& name)
    : G4VSensitiveDetector(name), fSimulationManager(simulationManager) {}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    // return value will always be ignored, its present for backwards compatibility (I guess)
    const auto volumeName = fSimulationManager->GetRestMetadata()
                                ->GetGeant4GeometryInfo()
                                .GetAlternativeNameFromGeant4PhysicalName(
                                    step->GetPreStepPoint()->GetPhysicalVolume()->GetName());

    const bool isGeantino = step->GetTrack()->GetParticleDefinition() == G4Geantino::Definition();

    if (isGeantino) {
        // Since geantinos don't deposit energy, the length traveled inside the volumes is stored as energy
        // (mm as keV)
        const auto length = step->GetStepLength() / CLHEP::mm;
        fSimulationManager->GetOutputManager()->AddSensitiveEnergy(length);
        return true;
    } else {
        auto energy = step->GetTotalEnergyDeposit() / keV;

        if (energy <= 0) {
            return true;
        }
        fSimulationManager->GetOutputManager()->AddSensitiveEnergy(energy);
        return true;
    }
}
