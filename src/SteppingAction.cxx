
#include "SteppingAction.h"

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4Track.h>

#include <G4DynamicParticle.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4SteppingManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

#include "SimulationManager.h"

using namespace std;

SteppingAction::SteppingAction(SimulationManager* simulationManager)
    : fSimulationManager(simulationManager) {}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* aStep) {
    TRestGeant4Track* restTrack = fSimulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = fSimulationManager->fRestGeant4Event;
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;

    const auto& geometryInfo = restG4Metadata->GetGeant4GeometryInfo();
    // Variables that describe a step are taken.
    const auto& volumeName = geometryInfo.GetAlternativeNameFromGeant4PhysicalName(
        (TString &&) aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());

    const auto& particle = aStep->GetTrack()->GetDefinition();
    const auto& particleID = particle->GetPDGEncoding();
    const auto& particleName = particle->GetParticleName();

    restG4Metadata->fGeant4PhysicsInfo.InsertParticleName(particleID, particleName);

    const auto process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    const auto& processID = process->GetProcessType() * 1000 + process->GetProcessSubType();
    const auto& processName = process->GetProcessName();

    restG4Metadata->fGeant4PhysicsInfo.InsertProcessName(processID, processName);

    const auto totalEnergyDeposit = aStep->GetTotalEnergyDeposit() / keV;
    const auto trackKineticEnergy = aStep->GetTrack()->GetKineticEnergy() / keV;

    auto sensitiveVolumeName =
        geometryInfo.GetAlternativeNameFromGeant4PhysicalName(restG4Metadata->GetSensitiveVolume());

    if (restTrack->GetParticleName() == "geantino" && sensitiveVolumeName.Data() == volumeName) {
        restG4Metadata->SetSaveAllEvents(true);
    }

    if (!aStep->GetPostStepPoint()->GetProcessDefinedStep()) {
        G4cout << endl;
        G4cout << endl;
        G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        G4cout << "An ERROR was detected on the G4Step process definition." << endl;
        G4cout << "This is a sign of problem in the restG4 particle definition" << endl;
        G4cout << endl;
        G4cout << "E.g. A definition of a gamma with 0keV energy" << endl;
        G4cout << endl;
        G4cout << "Please, review your TRestGeant4Metadata RML definition" << endl;
        G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        G4cout << endl;
        G4cout << endl;

        exit(0);
    }

    G4Track* aTrack = aStep->GetTrack();

    Double_t x = aTrack->GetPosition().x() / mm;
    Double_t y = aTrack->GetPosition().y() / mm;
    Double_t z = aTrack->GetPosition().z() / mm;

    if (restG4Metadata->GetSensitiveVolume() == volumeName) {
        restG4Event->AddEnergyToSensitiveVolume(totalEnergyDeposit);
    }

    const TVector3 hitPosition(x, y, z);
    const Double_t hitGlobalTime = aStep->GetPreStepPoint()->GetGlobalTime() / second;
    const G4ThreeVector& momentum = aStep->GetPreStepPoint()->GetMomentumDirection();
    const TVector3 momentumDirection = TVector3(momentum.x(), momentum.y(), momentum.z());  //.Unit();

    Int_t volume = -1;
    Bool_t alreadyStored = false;
    // We check if the hit must be stored and keep it on restG4Track
    for (int volID = 0; volID < restG4Metadata->GetNumberOfActiveVolumes(); volID++) {
        if (restG4Event->isVolumeStored(volID)) {
            if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) {
                G4cout << "Step volume :" << volumeName << "::("
                       << (G4String)restG4Metadata->GetActiveVolumeName(volID) << ")" << G4endl;
            }

            // We store the hit if we have activated in the config
            const Bool_t isActiveVolume = (volumeName == restG4Metadata->GetActiveVolumeName(volID));
            if (isActiveVolume) {
                volume = volID;
                if (restG4Metadata->GetVerboseLevel() >=
                    TRestStringOutput::REST_Verbose_Level::REST_Extreme) {
                    G4cout << "Storing hit" << G4endl;
                }
                restTrack->AddG4Hit(hitPosition, totalEnergyDeposit, hitGlobalTime, processID, volID,
                                    trackKineticEnergy, momentumDirection);
                alreadyStored = true;
                restG4Metadata->fGeant4GeometryInfo.InsertVolumeName(volID, volumeName);
            }
        }
    }

    // See issue #65.
    // If the radioactive decay occurs in a non-active volume then the id will be -1
    const Bool_t isDecay = (processName == (G4String) "RadioactiveDecay");
    if (!alreadyStored && isDecay) {
        restTrack->AddG4Hit(hitPosition, totalEnergyDeposit, hitGlobalTime, processID, volume,
                            trackKineticEnergy, momentumDirection);
    }
}