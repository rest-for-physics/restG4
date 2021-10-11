
#ifndef RunAction_h
#define RunAction_h 1

#include <G4UserRunAction.hh>
#include <globals.hh>
#include <iostream>

using namespace std;

class G4Run;
class PrimaryGeneratorAction;
class TRestGeant4Metadata;
class OutputManager;

class RunAction : public G4UserRunAction {
   public:
    RunAction(PrimaryGeneratorAction*);
    ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

   private:
    PrimaryGeneratorAction* fPrimary;
    TRestGeant4Metadata* fRestGeant4Metadata;
    OutputManager* fOutputManager;
};

#endif
