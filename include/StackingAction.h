
#ifndef REST_STACKINGACTION_H
#define REST_STACKINGACTION_H

#include <G4ParticleDefinition.hh>
#include <G4ThreeVector.hh>
#include <G4UserStackingAction.hh>
#include <globals.hh>
#include <map>
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

    struct SubEventPrimary {
        G4String name;
        G4double energy;
        G4ThreeVector position;
        G4ThreeVector direction;
    };

    const SubEventPrimary* GetSubEventPrimary(G4int parentID) const {
        const auto primary = fSubEventPrimaries.find(parentID);
        return primary == fSubEventPrimaries.end() ? nullptr : &primary->second;
    }

    inline std::set<const G4ParticleDefinition*> GetParticlesToIgnore() const { return fParticlesToIgnore; }
    inline void AddParticleToIgnore(const G4ParticleDefinition* particle) {
        fParticlesToIgnore.insert(particle);
    }

   private:
    SimulationManager* fSimulationManager;
    double fMaxAllowedLifetime;
    std::map<G4int, SubEventPrimary> fSubEventPrimaries;
    std::set<G4int> fKillDecaysFromParentID;
    G4String fMaxAllowedLifetimeWithUnit;

    std::set<const G4ParticleDefinition*> fParticlesToIgnore;
};

#endif  // REST_STACKINGACTION_H
