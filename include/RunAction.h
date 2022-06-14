
#ifndef RunAction_h
#define RunAction_h 1

#include <G4UserRunAction.hh>
#include <fstream>
#include <globals.hh>
#include <iostream>
#include <map>

class G4Run;
class PrimaryGeneratorAction;

class RunAction : public G4UserRunAction {
   public:
    RunAction();
    ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

   private:
    std::map<G4String, G4int> fParticleCount;
    std::map<G4String, G4double> fEmean;
    std::map<G4String, G4double> fEmin;
    std::map<G4String, G4double> fEmax;
};

#endif
