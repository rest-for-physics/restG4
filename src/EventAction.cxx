
#include "EventAction.h"

#include <OutputManager.h>
#include <TRestRun.h>
#include <spdlog/spdlog.h>

#include <G4Event.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4Threading.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>

#include "GlobalManager.h"

using namespace std;

EventAction::EventAction()
    : G4UserEventAction(),
      fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()),
      fOutputManager(OutputManager::Instance()) {
    fRestGeant4Metadata->isFullChainActivated();
}

EventAction::~EventAction() = default;

void EventAction::BeginOfEventAction(const G4Event* event) {
    spdlog::debug("EventAction::BeginOfEventAction <--- Begin of event {}", event->GetEventID());

    fOutputManager->UpdateEvent();

    auto eventID = event->GetEventID();
    int s = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() / 100;
    if ((s > 0 && (eventID + 1) % s == 0) || eventID == 0) {
        spdlog::info(
            "EventAction::BeginOfEventAction - RunID: {} ---> Begin of event {} ({:03.2f}%) - Number of "
            "entries saved this run: {}",
            G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID(), eventID,
            100 * float(eventID + 1) /
                static_cast<float>(G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()),
            GlobalManager::Instance()->GetEntries());
    }

    return;
}

void EventAction::EndOfEventAction(const G4Event* event) {
    // G4String energyWithUnits = G4BestUnit(fOutputManager->GetSensitiveVolumeEnergy() * CLHEP::keV,
    // "Energy");

    spdlog::debug("EventAction::EndOfEventAction <--- End of event {}", event->GetEventID());

    auto restRun = GlobalManager::Instance()->GetRestRun();

    fOutputManager->FinishAndSubmitEvent();

    if (!G4Threading::IsMultithreadedApplication() ||  //
        (G4Threading::IsMultithreadedApplication() && G4Threading::G4GetThreadId() == 0)) {
        GlobalManager::Instance()->FillEvents();
    }

    return;
}
