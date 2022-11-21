
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
#ifdef GEANT4_WITHOUT_G4RunManagerFactory  // For old Geant4 where the thread print gives segfault
    else {
        const int numberOfEventsToBePercent =
            G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() / 100;
        if ((restG4Metadata->PrintProgress() &&
             restG4Metadata->GetVerboseLevel() != TRestStringOutput::REST_Verbose_Level::REST_Silent) &&
            // print roughly every 1% of events or whenever 10 seconds without printing have elapsed
            (numberOfEventsToBePercent > 0 && (eventID + 1) % numberOfEventsToBePercent == 0)) {
            fSimulationManager->SyncStatsFromChild();
            G4cout << double(fSimulationManager->GetNumberOfProcessedEvents()) /
                          double(restG4Metadata->GetNumberOfEvents()) * 100
                   << "% - " << fSimulationManager->GetNumberOfProcessedEvents()
                   << " Events processed out of " << restG4Metadata->GetNumberOfEvents()
                   << " requested events ("
                   << fSimulationManager->GetNumberOfProcessedEvents() / fSimulationManager->GetElapsedTime()
                   << " per second). " << fSimulationManager->GetNumberOfStoredEvents() << " events stored ("
                   << fSimulationManager->GetNumberOfStoredEvents() / fSimulationManager->GetElapsedTime()
                   << " per second). " << fSimulationManager->GetElapsedTime() << " seconds elapsed" << endl;
        }
    }
#endif
}

void EventAction::EndOfEventAction(const G4Event*) {
    fSimulationManager->GetOutputManager()->FinishAndSubmitEvent();
}
