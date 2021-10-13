//
// Created by lobis on 12/10/2021.
//

#include <SteppingAction.h>
#include <TRestGeant4DataEvent.h>
#include <TRestGeant4DataSteps.h>
#include <TRestGeant4DataTrack.h>

#include <G4Event.hh>
#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include <G4ParticleDefinition.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4Scintillation.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4Transportation.hh>
#include <G4UnitsTable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VProcess.hh>
#include <G4VUserPrimaryVertexInformation.hh>

#include "GlobalManager.h"

using namespace std;

G4String GetVolumeName(G4VPhysicalVolume* volume) {
    if (!volume) {
        const G4String invalidVolumeName = "?";
        return invalidVolumeName;
    }
    auto defaultName = volume->GetName();

    /*
    auto lookupName = GlobalManager::Instance()->GetVolumeFromLookupTable(defaultName);
    if (lookupName != "") {
        return lookupName;
    }
    */
    return defaultName;
}

void TRestGeant4DataSteps::InsertStep(const G4Step* step) {
    const G4Track* track = step->GetTrack();

    const auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4String processName = "Init";
    G4String processTypeName = "Init";
    G4String targetNucleusName = "";
    G4int processID = -1;
    if (track->GetCurrentStepNumber() != 0) {
        // 0 = Init step (G4SteppingVerbose) process is not defined for this step
        processName = process->GetProcessName();
        processTypeName = G4VProcess::GetProcessTypeName(process->GetProcessType());
        processID = process->GetProcessType() * 1000 + process->GetProcessSubType();
        //
        if (process->GetProcessType() == G4ProcessType::fHadronic) {
            auto hadronicProcess = dynamic_cast<const G4HadronicProcess*>(process);
            G4Nucleus nucleus = *(hadronicProcess->GetTargetNucleus());
            auto isotope = nucleus.GetIsotope();
            if (isotope) {
                targetNucleusName = nucleus.GetIsotope()->GetName();
            }
        }
    }

    const auto volumeNamePre = GetVolumeName(step->GetPreStepPoint()->GetPhysicalVolume());
    auto volumeNamePost = (step->GetPostStepPoint()->GetPhysicalVolume()
                               ? GetVolumeName(step->GetPostStepPoint()->GetPhysicalVolume())
                               : "?");
    if (volumeNamePre == volumeNamePost) {
        volumeNamePost = "";
    }
    const auto volumeID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
    const auto& volumeName = volumeNamePre;

    fN += 1;

    fStepID.emplace_back(track->GetCurrentStepNumber());
    fPosition.emplace_back(track->GetPosition().x() / CLHEP::mm, track->GetPosition().y() / CLHEP::mm,
                           track->GetPosition().z() / CLHEP::mm);
    fMomentumDirection.emplace_back(track->GetMomentum().x() / CLHEP::keV,
                                    track->GetMomentum().y() / CLHEP::keV,
                                    track->GetMomentum().z() / CLHEP::keV);
    fTimeGlobal.emplace_back(track->GetGlobalTime() / CLHEP::ns);
    fEnergy.emplace_back(step->GetTotalEnergyDeposit() / CLHEP::keV);
    fKineticEnergy.emplace_back(step->GetPreStepPoint()->GetKineticEnergy() / CLHEP::keV);
    fKineticEnergy.emplace_back(step->GetPostStepPoint()->GetKineticEnergy() / CLHEP::keV);
    fLength.emplace_back(step->GetStepLength() / CLHEP::mm);
    fVolumeName.emplace_back(volumeName);
    fVolumeNamePost.emplace_back(volumeNamePost);
    fVolumeID.emplace_back(volumeID);
    fProcessName.emplace_back(processName);
    fProcessID.emplace_back(processID);
    fProcessType.emplace_back(processTypeName);
    fTargetNucleus.emplace_back(targetNucleusName);
    //
    G4String energyWithUnits = G4BestUnit(fEnergy.back() * CLHEP::keV, "Energy");
    /*
    spdlog::debug(
        "DataModelSteps::InsertStep - Step ID {} - process {} - energy deposited {} - volume {}{} -
    position "
        "(mm) ({:03.2f}, {:03.2f}, {:03.2f})",                                     //
        fStepID.back(), fProcessName.back(), energyWithUnits, fVolumeName.back(),  //
        (fVolumeNamePost.back().IsNull() ? "" : "->" + fVolumeNamePost.back()),    //
        fPosition.back().x(), fPosition.back().y(), fPosition.back().z()           //
    );
     */
}
