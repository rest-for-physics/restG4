
#ifndef REST_STACKINGACTION_H
#define REST_STACKINGACTION_H

#include <G4ParticleDefinition.hh>
#include <G4UserStackingAction.hh>
#include <globals.hh>
#include <set>

#include "SimulationManager.h"

class OutputManager;
class SimulationManager;

class StackingAction : public G4UserStackingAction {
   public:
    explicit StackingAction(SimulationManager*);
    ~StackingAction();

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    void NewStage();
    void PrepareNewEvent();

    inline std::set<const G4ParticleDefinition*> GetParticlesToIgnore() const { return fParticlesToIgnore; }
    inline void AddParticleToIgnore(const G4ParticleDefinition* particle) {
        fParticlesToIgnore.insert(particle);
    }

   private:
    SimulationManager* fSimulationManager;
    double fMaxAllowedLifetime;
    std::set<G4int> fNewSubEventFromParentID;
    std::set<G4int> fKillDecaysFromParentID;
    G4String fMaxAllowedLifetimeWithUnit;

    std::set<const G4ParticleDefinition*> fParticlesToIgnore;
};

#endif  // REST_STACKINGACTION_H
