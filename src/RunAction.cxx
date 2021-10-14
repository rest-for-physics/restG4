
#include "RunAction.h"

#include <TRestGeant4Metadata.h>
#include <spdlog/spdlog.h>

#include <G4PhysicalConstants.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>

#include "GlobalManager.h"
#include "OutputManager.h"
#include "PrimaryGeneratorAction.h"
#include "SteppingVerbose.h"

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

void RunAction::BeginOfRunAction(const G4Run* run) {
    spdlog::info("RunAction::BeginOfRunAction ---> Begin of run {}", run->GetRunID());

    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    auto steppingVerbose = ((SteppingVerbose*)G4VSteppingVerbose::GetInstance());
    // steppingVerbose->SetSteppingVerbose(1);
}

void RunAction::EndOfRunAction(const G4Run* run) {
    spdlog::info("RunAction::EndOfRunAction <--- End of run {} (thread {})", run->GetRunID(),
                 G4Threading::G4GetThreadId());

    if (G4Threading::IsMasterThread()) {
        spdlog::info("RunAction::EndOfRunAction <--- Writing events to file");
        GlobalManager::Instance()->WriteEvents();
    }
}
