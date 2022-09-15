
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
        G4cout << "========================== Begin of Run Action ========================" << endl;
        G4cout << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << " events to be simulated"
               << endl;
        G4cout << "=======================================================================" << endl;
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
