
#include "EventAction.h"

#include <TRestRun.h>

#include <G4Event.hh>
#include <Randomize.hh>
#include <fstream>
#include <iomanip>

extern TRestRun* restRun;
extern TRestGeant4Metadata* restG4Metadata;
extern TRestGeant4Event* restG4Event;
extern TRestGeant4Event* subRestG4Event;
extern TRestGeant4Track* restTrack;

using namespace std;

EventAction::EventAction() : G4UserEventAction() { restG4Metadata->isFullChainActivated(); }

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event* geant4_event) {
    G4int event_number = geant4_event->GetEventID();

    restG4Metadata->GetVerboseLevel();

    if (restG4Metadata->GetVerboseLevel() >= REST_Debug) {
        cout << "DEBUG: Start of event ID " << event_number << " (" << event_number + 1 << " of "
             << restG4Metadata->GetNumberOfEvents() << "). " << restRun->GetEntries() << " Events stored."
             << endl;
    } else if ((restG4Metadata->PrintProgress() || restG4Metadata->GetVerboseLevel() >= REST_Info) &&
               geant4_event->GetEventID() % 10000 == 0) {
        cout << "INFO: Start of event ID " << event_number << " (" << event_number + 1 << " of "
             << restG4Metadata->GetNumberOfEvents() << "). " << restRun->GetEntries() << " Events stored."
             << endl
             << endl;
    }

    restTrack->Initialize();

    restG4Event->SetID(event_number);
    restG4Event->SetOK(true);
    time_t system_time = time(nullptr);

    restG4Event->SetTime((Double_t)system_time);

    // Defining if the hits in a given volume will be stored
    for (int i = 0; i < restG4Metadata->GetNumberOfActiveVolumes(); i++) {
        Double_t rndNumber = G4UniformRand();

        if (restG4Metadata->GetStorageChance(i) >= rndNumber)
            restG4Event->ActivateVolumeForStorage(i);
        else
            restG4Event->DisableVolumeForStorage(i);
    }
}

void EventAction::EndOfEventAction(const G4Event* geant4_event) {
    G4int event_number = geant4_event->GetEventID();

    if (restG4Metadata->GetVerboseLevel() >= REST_Extreme) {
        restG4Event->PrintEvent();
    }

    Double_t minimum_energy_stored = restG4Metadata->GetMinimumEnergyStored();
    Double_t maximum_energy_stored = restG4Metadata->GetMaximumEnergyStored();

    int NSubs = SetTrackSubEventIDs();

    Bool_t is_sensitive = false;

    for (int subId = 0; subId < NSubs + 1; subId++) {
        FillSubEvent(subId);

        Double_t total_deposited_energy = subRestG4Event->GetTotalDepositedEnergy();
        Double_t sensitive_volume_deposited_energy = subRestG4Event->GetSensitiveVolumeEnergy();

        if (minimum_energy_stored < 0) minimum_energy_stored = 0;
        if (maximum_energy_stored == 0) maximum_energy_stored = total_deposited_energy + 1.;

        is_sensitive =
            (sensitive_volume_deposited_energy > 0 && total_deposited_energy > minimum_energy_stored &&
             total_deposited_energy < maximum_energy_stored) ||
            restG4Metadata->GetSaveAllEvents();

        if (restG4Metadata->GetVerboseLevel() >= REST_Info) {
            string debug_level = "INFO";
            if (restG4Metadata->GetVerboseLevel() >= REST_Debug) {
                debug_level = "DEBUG";
            }

            if (is_sensitive || restG4Metadata->GetVerboseLevel() >= REST_Debug) {
                cout << debug_level
                     << ": Energy deposited in ACTIVE and SENSITIVE volumes: " << total_deposited_energy
                     << " keV" << endl;
                cout << debug_level
                     << ": Energy deposited in SENSITIVE volume: " << sensitive_volume_deposited_energy
                     << " keV" << endl;
            }
        }

        if (is_sensitive) {
            sensitive_volume_hits_count += 1;

            // call `ReOrderTrackIds` which before was integrated into `FillSubEvent`
            // it takes a while to run so we only do it if we are going to save the event
            ReOrderTrackIds(subId);

            // fill analysis tree
            TRestAnalysisTree* analysis_tree = restRun->GetAnalysisTree();
            if (analysis_tree != nullptr) {
                analysis_tree->SetEventInfo(subRestG4Event);
                analysis_tree->Fill();
            } else {
                // analysis tree is not found (nullptr)
                if (restG4Metadata->GetVerboseLevel() >= REST_Warning) {
                    cout << "WARNING: Analysis tree is not found ('nullptr'). Cannot write event info"
                         << endl;
                }
            }
            // fill event tree
            TTree* event_tree = restRun->GetEventTree();
            if (event_tree != nullptr) {
                event_tree->Fill();
            } else {
                // event tree is not found (nullptr)
                if (restG4Metadata->GetVerboseLevel() >= REST_Warning) {
                    cout << "WARNING: Event tree is not found ('nullptr'). Cannot write event info" << endl;
                }
            }
        }
    }

    if (restG4Metadata->GetVerboseLevel() >= REST_Info) {
        string debug_level = "INFO";
        if (restG4Metadata->GetVerboseLevel() >= REST_Debug) {
            debug_level = "DEBUG";
        }

        if (is_sensitive || restG4Metadata->GetVerboseLevel() >= REST_Debug) {
            cout << debug_level
                 << ": Events depositing energy in sensitive volume: " << sensitive_volume_hits_count << "/"
                 << event_number + 1 << endl;
            cout << debug_level << ": End of event ID " << event_number << " (" << event_number + 1 << " of "
                 << restG4Metadata->GetNumberOfEvents() << ")" << endl;
            cout << endl;
        }
    }
}

/*
 * TODO:
 * this method takes a very long time and is acting as a bottleneck on some instances (high track
 * number), see if it can be optimised.
 * */
void EventAction::FillSubEvent(Int_t subId) {
    subRestG4Event->Initialize();
    subRestG4Event->ClearVolumes();

    subRestG4Event->SetID(restG4Event->GetID());
    subRestG4Event->SetSubID(subId);

    subRestG4Event->SetRunOrigin(restRun->GetRunNumber());
    subRestG4Event->SetSubRunOrigin(0);

    time_t systime = time(NULL);
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
        const auto& tck = restG4Event->GetTrack(n);

        if (tck.GetSubEventID() == subId)
            if (tck.GetNumberOfHits() > 0 || restG4Metadata->RegisterEmptyTracks()) {
                subRestG4Event->AddTrack(tck);
            }
    }

    if (restG4Metadata->isVolumeStored(restG4Metadata->GetSensitiveVolume())) {
        Int_t sensVolID = restG4Metadata->GetActiveVolumeID(restG4Metadata->GetSensitiveVolume());

        subRestG4Event->SetSensitiveVolumeEnergy(subRestG4Event->GetEnergyDepositedInVolume(sensVolID));
    }
}

void EventAction::ReOrderTrackIds(Int_t subId) {
    // We define as event timestamp the system time.
    // We will be always able to extract the global simulation time from Geant4 tracks.
    time_t systime = time(NULL);
    subRestG4Event->SetTimeStamp(systime);

    if (subId > 0) {
        for (int n = 0; n < restG4Event->GetNumberOfTracks(); n++) {
            const auto& tck = restG4Event->GetTrack(n);

            if (tck.GetSubEventID() == subId - 1)
                if (tck.isRadiactiveDecay()) subRestG4Event->SetSubEventTag(tck.GetParticleName());
        }
    }

    // Re-ordering track IDs
    Int_t lowestID = subRestG4Event->GetLowestTrackID();
    Int_t nTracks = subRestG4Event->GetNumberOfTracks();

    for (int i = 0; i < nTracks; i++) {
        const auto& tr = subRestG4Event->GetTrackPointer(i);
        tr->SetTrackID(tr->GetTrackID() - lowestID + 1);
        tr->SetParentID(tr->GetParentID() - lowestID + 1);
        if (tr->GetParentID() < 0) tr->SetParentID(0);
    }

    for (int i = 0; i < nTracks; i++) {
        TRestGeant4Track* tr = subRestG4Event->GetTrackPointer(i);
        Int_t id = tr->GetTrackID();

        if (id - i != 1) {
            // Changing track ids
            tr->SetTrackID(i + 1);
            for (int t = i + 1; t < subRestG4Event->GetNumberOfTracks(); t++) {
                TRestGeant4Track* tr2 = subRestG4Event->GetTrackPointer(t);
                if (tr2->GetTrackID() == i + 1) tr2->SetTrackID(id);
            }

            // Changing parent ids
            for (int t = 0; t < subRestG4Event->GetNumberOfTracks(); t++) {
                TRestGeant4Track* tr2 = subRestG4Event->GetTrackPointer(t);
                if (tr2->GetParentID() == id)
                    tr2->SetParentID(i + 1);
                else if (tr2->GetParentID() == i + 1)
                    tr2->SetParentID(id);
            }
        }
    }
}

int EventAction::SetTrackSubEventIDs() {
    Int_t nTracks = restG4Event->GetNumberOfTracks();
    Double_t timeDelay = restG4Metadata->GetSubEventTimeDelay();  // in unit us

    // reorder tracks
    std::map<int, TRestGeant4Track*> tracks;
    for (int n = 0; n < nTracks; n++) {
        TRestGeant4Track* track = restG4Event->GetTrackPointer(n);
        tracks[track->GetTrackID()] = track;
    }

    // scan backwards to the parent tracks and set the track's sub event id
    int max_subid = 0;
    for (auto iter = tracks.begin(); iter != tracks.end(); iter++) {
        TRestGeant4Track* track = iter->second;

        if (track->GetParentID() == 0) {
            track->SetSubEventID(0);
        } else {
            double tadd = 0;
            int parentid = track->GetParentID();
            TRestGeant4Track* ptrack = tracks[parentid];
            while (1) {
                if (ptrack != NULL) {
                    tadd += ptrack->GetTrackTimeLength();
                    if (tadd > timeDelay) {
                        int subid = ptrack->GetSubEventID() + 1;
                        track->SetSubEventID(subid);
                        if (max_subid < subid) {
                            max_subid = subid;
                        }
                        break;
                    } else {
                        if (ptrack->GetParentID() == 0) {
                            track->SetSubEventID(0);
                            break;
                        } else {
                            ptrack = tracks[ptrack->GetParentID()];
                            continue;
                        }
                    }
                } else {
                    cout << "error! parent track is null" << endl;
                    abort();
                }
            }
        }
    }

    return max_subid;
}
