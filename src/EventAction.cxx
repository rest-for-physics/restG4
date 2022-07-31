
#include "EventAction.h"

#include <TRestRun.h>

#include <G4Event.hh>
#include <G4RunManager.hh>
#include <Randomize.hh>
#include <fstream>
#include <iomanip>

#include "SimulationManager.h"

using namespace std;

EventAction::EventAction(SimulationManager* simulationManager)
    : G4UserEventAction(), fSimulationManager(simulationManager) {}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fSimulationManager->GetOutputManager()->BeginOfEventAction();

    const auto eventID = event->GetEventID();
    TRestRun* restRun = fSimulationManager->GetRestRun();

    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        G4cout << "DEBUG: Start of event ID " << eventID << " (" << eventID + 1 << " of "
               << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << "). "
               << restRun->GetEntries() << " Events stored" << endl;
    }
}

void EventAction::EndOfEventAction(const G4Event*) {
    fSimulationManager->GetOutputManager()->FinishAndSubmitEvent();
}
