
#include "RunAction.h"

#include <TRestRun.h>

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <iomanip>

#include "PrimaryGeneratorAction.h"
#include "SimulationManager.h"
#include "SteppingVerbose.h"

using namespace std;

RunAction::RunAction(SimulationManager* simulationManager)
    : G4UserRunAction(), fSimulationManager(simulationManager) {}

RunAction::~RunAction() = default;

void RunAction::BeginOfRunAction(const G4Run*) {
    if (G4Threading::IsMasterThread() || !G4Threading::IsMultithreadedApplication()) {
        G4cout << "========================== Begin of Run Action ========================" << G4endl;
        G4cout << "Events to be simulated: "
               << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << G4endl;
        if (fSimulationManager->GetRestMetadata()->GetNumberOfRequestedEntries() > 0) {
            G4cout << "Requested number of entries in file: "
                   << fSimulationManager->GetRestMetadata()->GetNumberOfRequestedEntries() << G4endl;
        }
        if (fSimulationManager->GetRestMetadata()->GetSimulationMaxTimeSeconds() > 0) {
            G4cout << "Maximum simulation time: "
                   << ToTimeStringLong(fSimulationManager->GetRestMetadata()->GetSimulationMaxTimeSeconds())
                   << G4endl;
        }
        G4cout << "=======================================================================" << G4endl;
    }

    fSimulationManager->BeginOfRunAction();

    auto steppingVerbose = ((SteppingVerbose*)G4VSteppingVerbose::GetInstance());
    steppingVerbose->SetSteppingVerbose(1);

    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void RunAction::EndOfRunAction(const G4Run*) {
    fSimulationManager->EndOfRunAction();

    TRestRun* restRun = fSimulationManager->GetRestRun();
    const auto metadata = fSimulationManager->GetRestMetadata();

    if (G4Threading::IsMasterThread() || !G4Threading::IsMultithreadedApplication()) {
        G4cout << "============================= Run Summary =============================" << endl;
        G4cout << restRun->GetEntries() << " events stored out of " << metadata->GetNumberOfEvents()
               << " simulated events" << endl;
        G4cout << "=======================================================================" << endl;
    }
}
