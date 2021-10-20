//
// Created by lobis on 12/10/2021.
//

#include <SteppingAction.h>
#include <TRestGeant4DataEvent.h>
#include <TRestGeant4DataSteps.h>
#include <TRestGeant4DataTrack.h>
#include <TRestRun.h>
#include <spdlog/spdlog.h>

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

// Constructors

TRestGeant4DataEvent::TRestGeant4DataEvent(const G4Event* event) : TRestGeant4DataEvent() {
    fRunOrigin = GlobalManager::Instance()->GetRestRun()->GetRunNumber();

    fRunID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    fEventID = event->GetEventID();
    fSubEventID = 0;

    for (int i = 0; i < event->GetNumberOfPrimaryVertex(); i++) {
        auto primaryVertex = event->GetPrimaryVertex();

        if (primaryVertex->GetNumberOfParticle() > 1) {
            cout << "multiple primary particles detected, but only recording first particle!" << endl;
            exit(1);
        }

        const auto& position = primaryVertex->GetPosition();
        fPrimaryPosition.emplace_back(position.x() / CLHEP::mm, position.y() / CLHEP::mm,
                                      position.z() / CLHEP::mm);
        auto primaryParticle = primaryVertex->GetPrimary();
        fPrimaryParticleName.emplace_back(primaryParticle->GetParticleDefinition()->GetParticleName());
        fPrimaryEnergy.emplace_back(primaryParticle->GetKineticEnergy() / CLHEP::keV);
        const auto& momentum = primaryParticle->GetMomentumDirection();
        fPrimaryDirection.emplace_back(momentum.x() / CLHEP::mm, momentum.y() / CLHEP::mm,
                                       momentum.z() / CLHEP::mm);
    }
}

bool IsValid(const G4Track* track) {
    auto excludedParticles = {"opticalphoton"};
    // optical photons take too much space to store them
    for (const auto& particleName : excludedParticles) {
        if (track->GetParticleDefinition()->GetParticleName() == particleName) {
            return false;
        }
    }
    return true;
}

bool IsValid(const G4Step* step) { return IsValid(step->GetTrack()); }

void TRestGeant4DataEvent::InsertStep(const G4Step* step) {
    if (!IsValid(step)) {
        return;
    }
    if (step->GetTrack()->GetCurrentStepNumber() == 0) {
        // initial step (from SteppingVerbose) is generated before TrackingAction can insert the first track
        fInitialStep = TRestGeant4DataSteps();
        fInitialStep.InsertStep(step);
    } else {
        fTracks.back().InsertStep(step);
    }
}

void TRestGeant4DataEvent::InsertTrack(const G4Track* track) {
    if (!IsValid(track)) {
        return;
    }

    if (fTracks.empty()) {
        // primary for the sub-event
        fSubEventPrimaryParticleName = track->GetParticleDefinition()->GetParticleName();
        fSubEventPrimaryEnergy = track->GetKineticEnergy() / CLHEP::keV;
        const auto& position = track->GetPosition();
        fSubEventPrimaryPosition =
            TVector3(position.x() / CLHEP::mm, position.y() / CLHEP::mm, position.z() / CLHEP::mm);
        const auto& momentum = track->GetMomentumDirection();
        fSubEventPrimaryDirection =
            TVector3(momentum.x() / CLHEP::mm, momentum.y() / CLHEP::mm, momentum.z() / CLHEP::mm);
    }
    fTracks.emplace_back(track);
    if (fInitialStep.GetNumberOfSteps() != 1) {
        spdlog::error("fInitialStep does not have exactly one step!");
        exit(1);
    }
    fTracks.back().UpdateSteps(fInitialStep);
}

void TRestGeant4DataEvent::UpdateTrack(const G4Track* track) {
    if (!IsValid(track)) {
        return;
    }
    fTracks.back().UpdateTrack(track);
}
