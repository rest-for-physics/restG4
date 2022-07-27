
#include "RunAction.h"

#include <PrimaryGeneratorAction.h>
#include <SteppingVerbose.h>
#include <TRestGeant4Metadata.h>
#include <TRestRun.h>

#include <G4PhysicalConstants.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <iomanip>

#include "SimulationManager.h"

using namespace std;

RunAction::RunAction(SimulationManager* simulationManager)
    : G4UserRunAction(), fSimulationManager(simulationManager) {}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*) {
    if (G4Threading::IsMasterThread() || !G4Threading::IsMultithreadedApplication()) {
        G4cout << "========================== Begin of Run Action ========================" << endl;
        G4cout << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << " events to be simulated"
               << endl;
        G4cout << "=======================================================================" << endl;
    }

    auto steppingVerbose = ((SteppingVerbose*)G4VSteppingVerbose::GetInstance());
    steppingVerbose->SetSteppingVerbose(1);

    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void RunAction::EndOfRunAction(const G4Run*) {
    fSimulationManager->EndOfRun();

    TRestRun* restRun = fSimulationManager->GetRestRun();
    const auto metadata = fSimulationManager->GetRestMetadata();

    if (G4Threading::IsMasterThread() || !G4Threading::IsMultithreadedApplication()) {
        G4cout << "============================= Run Summary =============================" << endl;
        G4cout << restRun->GetEntries() << " events stored out of " << metadata->GetNumberOfEvents()
               << " simulated events" << endl;
        G4cout << "=======================================================================" << endl;
    }

    // Sanity check
    if (metadata->GetNumberOfDesiredEntries() == 0 && metadata->GetSimulationMaxTimeSeconds() == 0 &&
        metadata->GetNumberOfEvents() != G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()) {
        G4cout << "FATAL ERROR: possible error when calculating number of processed events, please check "
                  "'RunAction::EndOfRunAction'"
               << endl;
        exit(1);
    }
}
