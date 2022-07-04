
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include <TH1D.h>
#include <TH2D.h>

#include <G4RunManager.hh>
#include <G4UserSteppingAction.hh>
#include <fstream>
#include <globals.hh>
#include <iostream>

class SimulationManager;

class SteppingAction : public G4UserSteppingAction {
   public:
    SteppingAction(SimulationManager*);
    ~SteppingAction();

    void UserSteppingAction(const G4Step*) override;

   private:
    SimulationManager* fSimulationManager;

    G4String nom_vol, nom_part, nom_proc;
    G4double ener_dep, eKin;

    G4ThreeVector previousDirection;
};
#endif
