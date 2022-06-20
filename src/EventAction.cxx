
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

    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;
    restG4Metadata->isFullChainActivated();
}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    const auto eventID = event->GetEventID();
    auto simulationManager = fSimulationManager;
    TRestRun* restRun = simulationManager->fRestRun;
    TRestGeant4Track* restTrack = simulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = simulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = simulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = simulationManager->fRestGeant4Metadata;

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        G4cout << "DEBUG: Start of event ID " << eventID << " (" << eventID + 1 << " of "
               << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << "). "
               << restRun->GetEntries() << " Events stored." << endl;
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

    restTrack->Initialize();


}

void EventAction::EndOfEventAction(const G4Event* event) {
    auto simulationManager = fSimulationManager;
    TRestRun* restRun = simulationManager->fRestRun;
    TRestGeant4Track* restTrack = simulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = simulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = simulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = simulationManager->fRestGeant4Metadata;

    G4int eventID = event->GetEventID();

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) {
        restG4Event->PrintEvent();
    }

    Double_t minimumEnergyStored = restG4Metadata->GetMinimumEnergyStored();
    Double_t maximumEnergyStored = restG4Metadata->GetMaximumEnergyStored();

    int nSubs = SetTrackSubEventIDs();

    Bool_t isSensitive = false;

    for (int subId = 0; subId < nSubs + 1; subId++) {
        FillSubEvent(subId);

        Double_t total_deposited_energy = subRestG4Event->GetTotalDepositedEnergy();
        Double_t sensitive_volume_deposited_energy = subRestG4Event->GetSensitiveVolumeEnergy();

        if (minimumEnergyStored < 0) minimumEnergyStored = 0;
        if (maximumEnergyStored == 0) maximumEnergyStored = total_deposited_energy + 1.;

        isSensitive =
            (sensitive_volume_deposited_energy > 0 && total_deposited_energy > minimumEnergyStored &&
             total_deposited_energy < maximumEnergyStored) ||
            restG4Metadata->GetSaveAllEvents();

        if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) {
            string debug_level = "INFO";
            if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
                debug_level = "DEBUG";
            }

            if (isSensitive ||
                restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
                G4cout << debug_level
                       << ": Energy deposited in ACTIVE and SENSITIVE volumes: " << total_deposited_energy
                       << " keV" << endl;
                G4cout << debug_level
                       << ": Energy deposited in SENSITIVE volume: " << sensitive_volume_deposited_energy
                       << " keV" << endl;
            }
        }

        if (isSensitive) {
            sensitive_volume_hits_count += 1;

            // call `ReOrderTrackIds` which before was integrated into `FillSubEvent`
            // it takes a while to run, so we only do it if we are going to save the event
            ReOrderTrackIds(subId);

            // fill analysis tree
            TRestAnalysisTree* analysis_tree = restRun->GetAnalysisTree();
            if (analysis_tree) {
                analysis_tree->SetEventInfo(subRestG4Event);
                analysis_tree->Fill();
            } else {
                // analysis tree is not found (nullptr)
                if (restG4Metadata->GetVerboseLevel() >=
                    TRestStringOutput::REST_Verbose_Level::REST_Warning) {
                    G4cout << "WARNING: Analysis tree is not found ('nullptr'). Cannot write event info"
                           << endl;
                }
            }
            // fill event tree
            TTree* eventTree = restRun->GetEventTree();
            if (eventTree) {
                eventTree->Fill();
            } else {
                // event tree is not found (nullptr)
                if (restG4Metadata->GetVerboseLevel() >=
                    TRestStringOutput::REST_Verbose_Level::REST_Warning) {
                    G4cout << "WARNING: Event tree is not found ('nullptr'). Cannot write event info" << endl;
                }
            }
        }
    }

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) {
        string debugLevel = "INFO";
        if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
            debugLevel = "DEBUG";
        }

        if (isSensitive ||
            restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
            G4cout << debugLevel
                   << ": Events depositing energy in sensitive volume: " << sensitive_volume_hits_count << "/"
                   << eventID + 1 << endl;
            G4cout << debugLevel << ": End of event ID " << eventID << " (" << eventID + 1 << " of "
                   << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << ")" << endl;
            G4cout << endl;
        }
    }
}

/*
 * TODO:
 * this method takes a very long time and is acting as a bottleneck on some instances (high track
 * number), see if it can be optimised.
 * */
void EventAction::FillSubEvent(Int_t subId) {
    TRestRun* restRun = fSimulationManager->fRestRun;
    TRestGeant4Track* restTrack = fSimulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = fSimulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = fSimulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;

    subRestG4Event->Initialize();
    subRestG4Event->ClearVolumes();

    subRestG4Event->SetID(restG4Event->GetID());
    subRestG4Event->SetSubID(subId);

    subRestG4Event->SetRunOrigin(restRun->GetRunNumber());
    subRestG4Event->SetSubRunOrigin(0);

    time_t systime = time(nullptr);
    subRestG4Event->SetTimeStamp((Double_t)systime);

    subRestG4Event->SetPrimaryEventOrigin(restG4Event->GetPrimaryEventOrigin());
    for (int n = 0; n < restG4Event->GetNumberOfPrimaries(); n++) {
        subRestG4Event->SetPrimaryEventParticleName(restG4Event->GetPrimaryEventParticleName(n));
        subRestG4Event->SetPrimaryEventDirection(restG4Event->GetPrimaryEventDirection(n));
        subRestG4Event->SetPrimaryEventEnergy(restG4Event->GetPrimaryEventEnergy(n));
    }

    for (int n = 0; n < restG4Event->GetNumberOfActiveVolumes(); n++) {
        subRestG4Event->AddActiveVolume((string)restG4Metadata->GetActiveVolumeName(n));
        if (restG4Event->isVolumeStored(n))
            subRestG4Event->ActivateVolumeForStorage(n);
        else
            subRestG4Event->DisableVolumeForStorage(n);
    }

    for (int n = 0; n < restG4Event->GetNumberOfTracks(); n++) {
        const auto& track = restG4Event->GetTrack(n);
        if (track.GetSubEventID() == subId) {
            if (track.GetNumberOfHits() > 0 || restG4Metadata->RegisterEmptyTracks()) {
                subRestG4Event->AddTrack(track);
            }
        }
    }

    if (restG4Metadata->isVolumeStored(restG4Metadata->GetSensitiveVolume())) {
        Int_t sensVolID = restG4Metadata->GetActiveVolumeID(restG4Metadata->GetSensitiveVolume());
        subRestG4Event->SetSensitiveVolumeEnergy(subRestG4Event->GetEnergyDepositedInVolume(sensVolID));
    }
}

void EventAction::ReOrderTrackIds(Int_t subId) {
    TRestRun* restRun = fSimulationManager->fRestRun;
    TRestGeant4Track* restTrack = fSimulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = fSimulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = fSimulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;

    // We define as event timestamp the system time.
    // We will be always able to extract the global simulation time from Geant4 tracks.
    time_t systime = time(nullptr);
    subRestG4Event->SetTimeStamp(systime);

    if (subId > 0) {
        for (int n = 0; n < restG4Event->GetNumberOfTracks(); n++) {
            const auto& track = restG4Event->GetTrack(n);
            if (track.GetSubEventID() == subId - 1) {
                if (track.isRadiactiveDecay()) {
                    subRestG4Event->SetSubEventTag(track.GetParticleName());
                }
            }
        }
    }

    // Re-ordering track IDs
    Int_t lowestID = subRestG4Event->GetLowestTrackID();
    Int_t nTracks = subRestG4Event->GetNumberOfTracks();

    for (int i = 0; i < nTracks; i++) {
        const auto& track = subRestG4Event->GetTrackPointer(i);
        track->SetTrackID(track->GetTrackID() - lowestID + 1);
        track->SetParentID(track->GetParentID() - lowestID + 1);
        if (track->GetParentID() < 0) {
            track->SetParentID(0);
        }
    }

    for (int i = 0; i < nTracks; i++) {
        TRestGeant4Track* track = subRestG4Event->GetTrackPointer(i);
        Int_t id = track->GetTrackID();

        if (id - i != 1) {
            // Changing track ids
            track->SetTrackID(i + 1);
            for (int t = i + 1; t < subRestG4Event->GetNumberOfTracks(); t++) {
                TRestGeant4Track* trackAux = subRestG4Event->GetTrackPointer(t);
                if (trackAux->GetTrackID() == i + 1) {
                    trackAux->SetTrackID(id);
                }
            }

            // Changing parent ids
            for (int t = 0; t < subRestG4Event->GetNumberOfTracks(); t++) {
                TRestGeant4Track* trackAux = subRestG4Event->GetTrackPointer(t);
                if (trackAux->GetParentID() == id) {
                    trackAux->SetParentID(i + 1);
                } else if (trackAux->GetParentID() == i + 1) {
                    trackAux->SetParentID(id);
                }
            }
        }
    }
}

int EventAction::SetTrackSubEventIDs() {
    TRestRun* restRun = fSimulationManager->fRestRun;
    TRestGeant4Track* restTrack = fSimulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = fSimulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = fSimulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->fRestGeant4Metadata;

    Int_t nTracks = restG4Event->GetNumberOfTracks();
    Double_t timeDelay = restG4Metadata->GetSubEventTimeDelay();  // in unit us

    // reorder tracks
    map<int, TRestGeant4Track*> tracks;
    for (int n = 0; n < nTracks; n++) {
        TRestGeant4Track* track = restG4Event->GetTrackPointer(n);
        tracks[track->GetTrackID()] = track;
    }

    // scan backwards to the parent tracks and set the track's sub event id
    int maxSubEventID = 0;
    for (auto iter = tracks.begin(); iter != tracks.end(); iter++) {
        TRestGeant4Track* track = iter->second;
        if (track->GetParentID() == 0) {
            track->SetSubEventID(0);
        } else {
            double trackTimeLengthTotal = 0;
            int parentID = track->GetParentID();
            TRestGeant4Track* trackAux = tracks[parentID];
            while (true) {
                if (trackAux) {
                    trackTimeLengthTotal += trackAux->GetTrackTimeLength();
                    if (trackTimeLengthTotal > timeDelay) {
                        int subEventIDAux = trackAux->GetSubEventID() + 1;
                        track->SetSubEventID(subEventIDAux);
                        if (maxSubEventID < subEventIDAux) {
                            maxSubEventID = subEventIDAux;
                        }
                        break;
                    } else {
                        if (trackAux->GetParentID() == 0) {
                            track->SetSubEventID(0);
                            break;
                        } else {
                            trackAux = tracks[trackAux->GetParentID()];
                            continue;
                        }
                    }
                } else {
                    G4cout << "error! parent track is null" << endl;
                    abort();
                }
            }
        }
    }

    return maxSubEventID;
}
