
#include "SteppingAction.h"

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsInfo.h>
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

SteppingAction::SteppingAction(SimulationManager* simulationManager) : fSimulationManager(simulationManager) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;
    Int_t& biasing = fSimulationManager->fBiasing;

    if (biasing > 1) restBiasingVolume = restG4Metadata->GetBiasingVolume(biasing - 1);
}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* aStep) {
    TRestRun* restRun = fSimulationManager->fRestRun;
    TRestGeant4Track* restTrack = fSimulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = fSimulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = fSimulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;
    TRestGeant4PhysicsLists* restPhysList = fSimulationManager->fRestGeant4PhysicsLists;
    Int_t& biasing = fSimulationManager->fBiasing;

    const auto& geometryInfo = restG4Metadata->GetGeant4GeometryInfo();
    // Variables that describe a step are taken.
    const auto& volumeName = geometryInfo.GetAlternativeNameFromGeant4PhysicalName(
        (TString &&) aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());

    const auto& particle = aStep->GetTrack()->GetDefinition();
    const auto& particleID = particle->GetPDGEncoding();
    const auto& particleName = particle->GetParticleName();

    restG4Metadata->fGeant4PhysicsInfo.InsertParticleName(particleID, particleName);

    const auto process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    const auto& processID = TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(process);
    const auto& processName = process->GetProcessName();

    restG4Metadata->fGeant4PhysicsInfo.InsertProcessName(processID, processName);

    ener_dep = aStep->GetTotalEnergyDeposit();
    eKin = aStep->GetTrack()->GetKineticEnergy() / keV;

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
    parentID = aTrack->GetParentID();
    trackID = aTrack->GetTrackID();

    Double_t x = aTrack->GetPosition().x() / mm;
    Double_t y = aTrack->GetPosition().y() / mm;
    Double_t z = aTrack->GetPosition().z() / mm;

    //   G4cout << "Step direction : " << aTrack->GetMomentumDirection() <<
    //   G4endl;

    if (biasing > 0) {
        // In biasing mode we do not store hits. Just check if we observe a gamma
        // inside the volume
        if (restBiasingVolume.isInside(x, y, z) && particleName == "gamma") {
            Double_t eKinetic = aStep->GetPreStepPoint()->GetKineticEnergy() / keV;

            // we add the gamma energy to the energy spectrum
            if (eKinetic > restBiasingVolume.GetMinEnergy() && eKinetic < restBiasingVolume.GetMaxEnergy()) {
                G4ThreeVector position = aStep->GetPreStepPoint()->GetPosition();
                G4ThreeVector positionNorm = -aStep->GetPreStepPoint()->GetPosition().unit();
                G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentumDirection();

                Double_t angle;
                if (restBiasingVolume.GetBiasingVolumeType() == "virtualBox") {
                    if (absDouble(position.x()) > absDouble(position.y()) &&
                        absDouble(position.x()) > absDouble(position.z())) {
                        // launching event from side x
                        if (position.x() > 0)
                            positionNorm.set(-1, 0, 0);
                        else
                            positionNorm.set(1, 0, 0);

                    } else if (absDouble(position.y()) > absDouble(position.x()) &&
                               absDouble(position.y()) > absDouble(position.z())) {
                        // launching event from side y
                        if (position.y() > 0)
                            positionNorm.set(0, -1, 0);
                        else
                            positionNorm.set(0, 1, 0);
                    } else {
                        // launching event from side y
                        if (position.z() > 0)
                            positionNorm.set(0, 0, -1);
                        else
                            positionNorm.set(0, 0, 1);
                    }
                }
                angle = momentum.angle(positionNorm);

                fBiasingSpectrum->Fill(eKinetic);
                fAngularDistribution->Fill(angle);

                if (restBiasingVolume.GetBiasingVolumeType() == "virtualSphere")
                    fSpatialDistribution->Fill(position.getPhi(), position.getTheta());
                else
                    fSpatialDistribution->Fill(position.x(), position.y());
            }

            // and we cancel the track
            //           aTrack->SetTrackStatus( fStopAndKill );
            aTrack->SetTrackStatus(fKillTrackAndSecondaries);

            // previousDirection = momentum;
        }

    } else {
        if (restG4Metadata->GetSensitiveVolume() == volumeName) {
            restG4Event->AddEnergyToSensitiveVolume(ener_dep / keV);
        }

        TVector3 hitPosition(x, y, z);
        Double_t hit_global_time = aStep->GetPreStepPoint()->GetGlobalTime() / second;
        G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentumDirection();
        TVector3 momentumDirection = TVector3(momentum.x(), momentum.y(), momentum.z());  //.Unit();

        Int_t volume = -1;
        Bool_t alreadyStored = false;
        // We check if the hit must be stored and keep it on restG4Track
        for (int volID = 0; volID < restG4Metadata->GetNumberOfActiveVolumes(); volID++) {
            if (restG4Event->isVolumeStored(volID)) {
                if (restG4Metadata->GetVerboseLevel() >=
                    TRestStringOutput::REST_Verbose_Level::REST_Extreme) {
                    G4cout << "Step volume :" << volumeName << "::("
                           << (G4String)restG4Metadata->GetActiveVolumeName(volID) << ")" << G4endl;
                }

                // We store the hit if we have activated in the config
                Bool_t isActiveVolume = (volumeName == restG4Metadata->GetActiveVolumeName(volID));

                if (isActiveVolume) {
                    volume = volID;
                    if (restG4Metadata->GetVerboseLevel() >=
                        TRestStringOutput::REST_Verbose_Level::REST_Extreme)
                        G4cout << "Storing hit" << G4endl;
                    restTrack->AddG4Hit(hitPosition, ener_dep / keV, hit_global_time, processID, volID, eKin,
                                        momentumDirection);
                    alreadyStored = true;
                    restG4Metadata->fGeant4GeometryInfo.InsertVolumeName(volID, volumeName);
                }
            }
        }

        // See issue #65.
        // If the radioactive decay occurs in a non-active volume then the id will be -1
        Bool_t isDecay = (processName == (G4String) "RadioactiveDecay");
        if (!alreadyStored && isDecay)
            restTrack->AddG4Hit(hitPosition, ener_dep / keV, hit_global_time, processID, volume, eKin,
                                momentumDirection);
    }
}

Int_t TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(const G4VProcess* process) {
    // This method should return a unique ID for each process
    return process->GetProcessType() * 1000 + process->GetProcessSubType();
}