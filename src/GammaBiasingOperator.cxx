//
// Created by lobis on 2/22/2023.
//

#include "GammaBiasingOperator.h"

#include "G4BiasingProcessInterface.hh"
#include "G4GenericMessenger.hh"
#include "GammaBiasingOperation.h"

GammaBiasingOperator::GammaBiasingOperator(int splittingFactor, bool biasOnlyOnce,
                                           const TVector3& biasingCenter)
    : G4VBiasingOperator("BremSplittingOperator"),
      fSplittingFactor(splittingFactor),
      fBiasOnlyOnce(biasOnlyOnce),
      fBiasingCenter(biasingCenter) {
    fBremSplittingOperation =
        new GammaBiasingOperation("BremSplittingOperation", splittingFactor, biasingCenter);
}

void GammaBiasingOperator::StartRun() { fBremSplittingOperation->SetSplittingFactor(fSplittingFactor); }

void GammaBiasingOperator::StartTracking(const G4Track*) { fNInteractions = 0; }

G4VBiasingOperation* GammaBiasingOperator::ProposeFinalStateBiasingOperation(
    const G4Track*, const G4BiasingProcessInterface*) {
    if (fBiasOnlyOnce && (fNInteractions > 0)) return nullptr;

    fNInteractions++;

    return fBremSplittingOperation;
}
