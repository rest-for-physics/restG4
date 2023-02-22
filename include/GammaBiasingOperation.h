//
// Created by lobis on 2/22/2023.
//

#ifndef REST_GAMMABIASINGOPERATION_H
#define REST_GAMMABIASINGOPERATION_H

#include <TVector3.h>

#include "G4BiasingProcessInterface.hh"
#include "G4ParticleChange.hh"
#include "G4VBiasingOperation.hh"

class GammaBiasingOperation : public G4VBiasingOperation {
   public:
    GammaBiasingOperation(const G4String& name, G4int splittingFactor, const TVector3& biasingCenter);

    virtual ~GammaBiasingOperation();

   public:
    virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw(
        const G4BiasingProcessInterface*, G4ForceCondition&) {
        return 0;
    }

    virtual G4VParticleChange* ApplyFinalStateBiasing(const G4BiasingProcessInterface*, const G4Track*,
                                                      const G4Step*, G4bool&);

    virtual G4double DistanceToApplyOperation(const G4Track*, G4double, G4ForceCondition*) { return DBL_MAX; }

    virtual G4VParticleChange* GenerateBiasingFinalState(const G4Track*, const G4Step*) { return 0; }

   public:
    void SetSplittingFactor(G4int splittingFactor) { fSplittingFactor = splittingFactor; }

    G4int GetSplittingFactor() const { return fSplittingFactor; }

   private:
    G4int fSplittingFactor;
    G4ParticleChange fParticleChange;
    TVector3 fBiasingCenter;
};

#endif  // REST_GAMMABIASINGOPERATION_H
