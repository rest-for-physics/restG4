
#include "RunAction.h"

#include <TRestGeant4Metadata.h>

#include <G4PhysicalConstants.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <iomanip>

#include "GlobalManager.h"
#include "OutputManager.h"
#include "PrimaryGeneratorAction.h"

RunAction::RunAction()
    : G4UserRunAction(),
      fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()),
      fOutputManager(nullptr) {
    fPrimary = (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

    if (!G4Threading::IsMasterThread() || !G4Threading::IsMultithreadedApplication()) {
        // RunAction is also instanced in the master thread but does not work so no need for output
        fOutputManager = OutputManager::Instance();
    }
}

RunAction::~RunAction() = default;

void RunAction::BeginOfRunAction(const G4Run*) {
    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void RunAction::EndOfRunAction(const G4Run* run) {
    if (G4Threading::IsMasterThread()) {
        // spdlog::info("RunAction::EndOfRunAction <--- Writing events to file");
        GlobalManager::Instance()->WriteEvents();
    }

    return;

    G4int nbEvents = run->GetNumberOfEvent();
    if (nbEvents == 0) {
        return;
    }

    G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4String partName = particle->GetParticleName();
    // G4double eprimary = fPrimary->GetParticleGun()->GetParticleEnergy();

    G4cout << "======================== run summary ======================";
    G4cout << "\n" << nbEvents << " Events simulated\n";
    G4cout << "===========================================================";
    G4cout << G4endl;
}
