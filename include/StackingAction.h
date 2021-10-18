//
// Created by lobis on 12/10/2021.
//

#ifndef REST_STACKINGACTION_H
#define REST_STACKINGACTION_H

#include <TRestGeant4Metadata.h>

#include <G4ParticleDefinition.hh>
#include <G4UserStackingAction.hh>
#include <globals.hh>
#include <set>

class OutputManager;

class StackingAction : public G4UserStackingAction {
   public:
    StackingAction();
    ~StackingAction();

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*) override;
    void NewStage() override;

    inline std::set<const G4ParticleDefinition*> GetParticlesToIgnore() const { return fParticlesToIgnore; }
    inline void AddParticleToIgnore(const G4ParticleDefinition* particle) {
        fParticlesToIgnore.insert(particle);
    }

   private:
    OutputManager* fOutputManager;
    TRestGeant4Metadata* fRestGeant4Metadata;

    G4double fMaxAllowedLifetime;
    G4String fMaxAllowedLifetimeWithUnit;

    std::set<const G4ParticleDefinition*> fParticlesToIgnore;

    Bool_t fFullChain = true;
};
#endif  // REST_STACKINGACTION_H
