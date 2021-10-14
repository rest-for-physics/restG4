//
// Created by lobis on 10/11/2021.
//

#ifndef REST_OUTPUTMANAGER_H
#define REST_OUTPUTMANAGER_H

#include <TFile.h>
#include <TRestGeant4DataEvent.h>
#include <TString.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "globals.hh"

class G4Event;
class G4Track;
class G4Step;

class OutputManager {
   public:
    static OutputManager* Instance();
    inline ~OutputManager() = default;

    void UpdateEvent();
    void FinishAndSubmitEvent();

    void RecordTrack(const G4Track*);
    void UpdateTrack(const G4Track*);
    void RecordStep(const G4Step*);

    inline Float_t GetSensitiveVolumeEnergy() const { return fSensitiveEnergyTotal; }
    inline void AddSensitiveEnergy(double energy) { fEvent->AddSensitiveVolumeEnergy(energy); }

   private:
    inline OutputManager() = default;
    static thread_local OutputManager* pinstance_;

    std::unique_ptr<TRestGeant4DataEvent> fEvent{};

    double fSensitiveEnergyTotal{};

    bool IsEmptyEvent() const;
    bool IsValidEvent() const;
};

#endif  // REST_OUTPUTMANAGER_H
