//
// Created by lobis on 10/11/2021.
//

#include "OutputManager.h"

#include <G4RunManager.hh>
#include <G4Threading.hh>
#include <G4Track.hh>
#include <G4UserEventAction.hh>

#include "GlobalManager.h"

using namespace std;

thread_local OutputManager* OutputManager::pinstance_ = nullptr;

OutputManager* OutputManager::Instance() {
    if (G4Threading::IsMasterThread() && G4Threading::IsMultithreadedApplication()) {
        cout << "OutputManager::Instance() - Thread local instance should never be invoked from master "
                "thread in "
                "a MT application"
             << endl;
    }
    if (pinstance_ == nullptr) {
        pinstance_ = new OutputManager();
    }
    return pinstance_;
}

void OutputManager::UpdateEvent() {
    auto event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    // fEvent = make_unique<TRestGeant4Event>(event);
}