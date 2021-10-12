
#include "SteppingAction.h"

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4Track.h>

#include <G4DynamicParticle.hh>
#include <G4OpticalPhoton.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4SteppingManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

#include "GlobalManager.h"
#include "OutputManager.h"

extern TRestGeant4Event* restG4Event;
extern TRestGeant4Track* restTrack;
extern Int_t biasing;

SteppingAction::SteppingAction()
    : fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()),
      fOutputManager(OutputManager::Instance()) {
    if (biasing > 1) restBiasingVolume = fRestGeant4Metadata->GetBiasingVolume(biasing - 1);
}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto track = step->GetTrack();

    // optical photons
    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        // TODO: Do something with optical photons
    }

    fOutputManager->RecordStep(step);

    return;

    // Variables that describe a step are taken.
    nom_vol = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    nom_part = step->GetTrack()->GetDefinition()->GetParticleName();

    ener_dep = step->GetTotalEnergyDeposit();
    eKin = step->GetTrack()->GetKineticEnergy() / keV;

    if (restTrack->GetParticleName() == "geantino" &&
        (G4String)fRestGeant4Metadata->GetSensitiveVolume() == nom_vol) {
        fRestGeant4Metadata->SetSaveAllEvents(true);
    }

    if (!step->GetPostStepPoint()->GetProcessDefinedStep()) {
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

    nom_proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

    G4Track* aTrack = step->GetTrack();
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
        if (restBiasingVolume.isInside(x, y, z) && nom_part == "gamma") {
            Double_t eKinetic = step->GetPreStepPoint()->GetKineticEnergy() / keV;

            // we add the gamma energy to the energy spectrum
            if (eKinetic > restBiasingVolume.GetMinEnergy() && eKinetic < restBiasingVolume.GetMaxEnergy()) {
                G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
                G4ThreeVector positionNorm = -step->GetPreStepPoint()->GetPosition().unit();
                G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentumDirection();

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
        if ((G4String)fRestGeant4Metadata->GetSensitiveVolume() == nom_vol) {
            restG4Event->AddEnergyToSensitiveVolume(ener_dep / keV);
        }

        TVector3 hitPosition(x, y, z);
        Int_t pcsID = restTrack->GetProcessID(nom_proc);
        Double_t hit_global_time = step->GetPreStepPoint()->GetGlobalTime() / second;
        G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentumDirection();
        TVector3 momentumDirection = TVector3(momentum.x(), momentum.y(), momentum.z());  //.Unit();

        Int_t volume = -1;
        Bool_t alreadyStored = false;
        // We check if the hit must be stored and keep it on restG4Track
        for (int volID = 0; volID < fRestGeant4Metadata->GetNumberOfActiveVolumes(); volID++) {
            if (restG4Event->isVolumeStored(volID)) {
                if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Extreme)
                    G4cout << "Step volume :" << nom_vol << "::("
                           << (G4String)fRestGeant4Metadata->GetActiveVolumeName(volID) << ")" << G4endl;

                // We store the hit if we have activated in the config
                Bool_t isActiveVolume =
                    (nom_vol == (G4String)fRestGeant4Metadata->GetActiveVolumeName(volID));

                if (isActiveVolume) {
                    volume = volID;
                    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Extreme)
                        G4cout << "Storing hit" << G4endl;
                    restTrack->AddG4Hit(hitPosition, ener_dep / keV, hit_global_time, pcsID, volID, eKin,
                                        momentumDirection);
                    alreadyStored = true;
                }
            }
        }

        // See issue #65.
        // If the radiactive decay occurs in a non active volume then the id will be -1
        Bool_t isDecay = (nom_proc == (G4String) "RadioactiveDecay");
        if (!alreadyStored && isDecay)
            restTrack->AddG4Hit(hitPosition, ener_dep / keV, hit_global_time, pcsID, volume, eKin,
                                momentumDirection);
    }
}

G4TrackVector* SteppingAction::GetfSecondary() {
    G4SteppingManager* steppingManager = fpSteppingManager;
    return steppingManager->GetfSecondary();
}
