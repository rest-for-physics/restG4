
#ifndef REST_ACTIONINITIALIZATION_H
#define REST_ACTIONINITIALIZATION_H

#include <G4VUserActionInitialization.hh>

class DetectorConstruction;
class G4VSteppingVerbose;
class SimulationManager;

class ActionInitialization : public G4VUserActionInitialization {
   public:
    ActionInitialization(SimulationManager*);
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

    // virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;

   private:
    SimulationManager* fSimulationManager;
};

#endif  // REST_ACTIONINITIALIZATION_H
