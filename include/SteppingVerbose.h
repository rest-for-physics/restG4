//
// Created by lobis on 6/20/2022.
//

#ifndef REST_STEPPINGVERBOSE_H
#define REST_STEPPINGVERBOSE_H

#include <G4SteppingManager.hh>
#include <G4SteppingVerbose.hh>
#include <G4VSteppingVerbose.hh>

class SimulationManager;
class G4Step;

class SteppingVerbose : public G4SteppingVerbose {
   public:
    SteppingVerbose(SimulationManager*);
    ~SteppingVerbose();

    virtual void TrackingStarted();
    virtual void StepInfo();

    int GetSteppingVerbose() { return fManager->GetverboseLevel(); }
    void SetSteppingVerbose(int level) { fManager->SetVerboseLevel(level); }

   private:
    SimulationManager* fSimulationManager;
};

#endif  // REST_STEPPINGVERBOSE_H
