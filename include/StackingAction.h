
#ifndef REST_STACKINGACTION_H
#define REST_STACKINGACTION_H

#include <G4ParticleDefinition.hh>
#include <G4UserStackingAction.hh>
#include <globals.hh>
#include <set>

class OutputManager;
class SimulationManager;

class StackingAction : public G4UserStackingAction {
   public:
    explicit StackingAction(SimulationManager*);
    ~StackingAction();

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    void NewStage();

    inline std::set<const G4ParticleDefinition*> GetParticlesToIgnore() const { return fParticlesToIgnore; }
    inline void AddParticleToIgnore(const G4ParticleDefinition* particle) {
        fParticlesToIgnore.insert(particle);
    }

   private:
    SimulationManager* fSimulationManager;

    const G4double fMaxAllowedLifetime;
    G4String fMaxAllowedLifetimeWithUnit;

    std::set<const G4ParticleDefinition*> fParticlesToIgnore;
};

#endif  // REST_STACKINGACTION_H
