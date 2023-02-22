//
// Created by lobis on 2/22/2023.
//

#include "GammaBiasingOperation.h"
#include "GammaBiasingOperator.h"

#include "G4BiasingProcessInterface.hh"
#include "G4GenericMessenger.hh"

GammaBiasingOperator::GammaBiasingOperator()
        : G4VBiasingOperator("BremSplittingOperator"),
          fSplittingFactor(1),
          fBiasOnlyOnce(true) {
    fBremSplittingOperation = new GammaBiasingOperation("BremSplittingOperation");
}


void GammaBiasingOperator::StartRun() {
    fBremSplittingOperation->SetSplittingFactor(fSplittingFactor);
}

void GammaBiasingOperator::StartTracking(const G4Track *) {
    fNInteractions = 0;
}

G4VBiasingOperation *
GammaBiasingOperator::ProposeFinalStateBiasingOperation(const G4Track *,
                                                        const G4BiasingProcessInterface *) {
    if (fBiasOnlyOnce && (fNInteractions > 0)) return nullptr;

    fNInteractions++;

    return fBremSplittingOperation;
}