
#include "SensitiveDetector.h"

#include <G4Geantino.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>

#include "SimulationManager.h"

using namespace std;

SensitiveDetector::SensitiveDetector(SimulationManager* simulationManager, const G4String& name)
    : G4VSensitiveDetector(name), fSimulationManager(simulationManager) {}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    // return value will always be ignored, its present for backwards compatibility (I guess)
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& geometryInfo = restG4Metadata->GetGeant4GeometryInfo();
    const auto volumeName = geometryInfo.GetAlternativeNameFromGeant4PhysicalName(
        step->GetPreStepPoint()->GetPhysicalVolume()->GetName());

    /* Sensitive detector are assigned to logical volumes, but several physical volumes
    can share the same logical volume. Now we get the physical volume where the pre-step
    is located and filter if it is one of the declared sensitive volumes by the user*/

    // Get the full name (path) of the physical volume which uniquely identifies it
    auto th = step->GetPreStepPoint()->GetTouchable();
    G4int depth = th->GetHistoryDepth();
    G4String geant4path = "";
    if (depth == 0) { // it is the world volume
        geant4path = th->GetVolume()->GetName();
    }
    for (G4int i = 1; i <= depth; ++i) {  // start from 1 to skip world volume
        // Move the touchable to level i (0 = current volume, depth = world)
        G4VPhysicalVolume* pv = th->GetVolume(depth - i);
        if (pv) {
            if (geant4path != "") {
                geant4path += geometryInfo.GetPathSeparator().Data();
            }
            geant4path += pv->GetName();
        }
    }
    // convert to the names used in gdml (due to assemblies)
    const auto physicalFullName = geometryInfo.GetAlternativePathFromGeant4Path(geant4path);

    // search and filter if this physical volume is sensitive
    auto sensitiveVolumes = restG4Metadata->GetSensitiveVolumes();
    if (std::find(sensitiveVolumes.begin(), sensitiveVolumes.end(), physicalFullName) ==
        sensitiveVolumes.end()) {
        return false;
    }

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
