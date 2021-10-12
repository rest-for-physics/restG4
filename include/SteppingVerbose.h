//
// Created by lobis on 12/10/2021.
//

#ifndef REST_STEPPINGVERBOSE_H
#define REST_STEPPINGVERBOSE_H

#include <G4SteppingManager.hh>
#include <G4SteppingVerbose.hh>
#include <G4VSteppingVerbose.hh>

class OutputManager;
class G4Step;

class SteppingVerbose : public G4SteppingVerbose {
   public:
    SteppingVerbose();
    ~SteppingVerbose();

    virtual void TrackingStarted();
    virtual void StepInfo();

    int GetSteppingVerbose() { return fManager->GetverboseLevel(); }
    void SetSteppingVerbose(int level) { fManager->SetVerboseLevel(level); }

   private:
    OutputManager* output;
};

#endif  // REST_STEPPINGVERBOSE_H
