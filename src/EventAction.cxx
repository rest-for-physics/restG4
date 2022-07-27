
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
    : G4UserEventAction(), fSimulationManager(simulationManager) {
    fTimer.Start();
}

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
    } else {
        const int numberOfEventsToBePercent =
            G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() / 100;
        if ((restG4Metadata->PrintProgress() ||
             restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Essential) &&
            // print roughly every 1% of events or whenever 10 seconds without printing have elapsed
            ((numberOfEventsToBePercent > 0 && (eventID + 1) % numberOfEventsToBePercent == 0) ||
             fTimer.RealTime() > 10.0)) {
            fTimer.Start();
            G4cout << "ESSENTIAL: Start of event ID " << eventID << " (" << eventID + 1 << " of "
                   << restG4Metadata->GetNumberOfEvents() << "). " << restRun->GetEntries()
                   << " Events stored" << endl;
        } else {
            fTimer.Continue();
        }
    }
}

void EventAction::EndOfEventAction(const G4Event*) {
    fSimulationManager->GetOutputManager()->FinishAndSubmitEvent();

    fSimulationManager->WriteEvents();
}
