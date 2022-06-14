
#include "RunAction.h"

#include <PrimaryGeneratorAction.h>
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
    G4cout << "========================== Begin of Run Action ========================" << endl;
    G4cout << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << " events to be simulated"
           << endl;
    G4cout << "=======================================================================" << endl;
    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void RunAction::EndOfRunAction(const G4Run*) {
    TRestRun* restRun = fSimulationManager->fRestRun;

    G4cout << "============================= Run Summary =============================" << endl;
    G4cout << restRun->GetEntries() << " events stored out of "
           << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << " simulated events" << endl;
    G4cout << "=======================================================================" << endl;
}
