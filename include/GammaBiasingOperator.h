//
// Created by lobis on 2/22/2023.
//

#ifndef REST_GAMMABIASINGOPERATOR_H
#define REST_GAMMABIASINGOPERATOR_H

#include <TVector3.h>

#include "G4VBiasingOperator.hh"

class GammaBiasingOperation;

class GammaBiasingOperator : public G4VBiasingOperator {
   public:
    GammaBiasingOperator(int splittingFactor, bool biasOnlyOnce, const TVector3& biasingCenter);

    virtual ~GammaBiasingOperator() {}

   public:
    virtual void StartRun();

    virtual void StartTracking(const G4Track* track);

    TVector3 GetBiasingCenter() const { return fBiasingCenter; }

   private:
    virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(const G4Track*,
                                                                   const G4BiasingProcessInterface*) {
        return 0;
    }

    virtual G4VBiasingOperation* ProposeOccurenceBiasingOperation(const G4Track*,
                                                                  const G4BiasingProcessInterface*) {
        return 0;
    }

    virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation(
        const G4Track* track, const G4BiasingProcessInterface* callingProcess);

   private:
    using G4VBiasingOperator::OperationApplied;

   private:
    GammaBiasingOperation* fBremSplittingOperation;
    G4int fSplittingFactor;
    G4bool fBiasOnlyOnce;
    G4int fNInteractions;
    TVector3 fBiasingCenter;
};

#endif  // REST_GAMMABIASINGOPERATOR_H
