//
// Created by lobis on 2/22/2023.
//

#include "GammaBiasingOperation.h"
#include "GammaBiasingOperator.h"

#include "G4BiasingProcessInterface.hh"
#include "G4GenericMessenger.hh"

GammaBiasingOperator::GammaBiasingOperator()
        : G4VBiasingOperator("BremSplittingOperator"),
          fSplittingFactor(10),
          fBiasPrimaryOnly(false),
          fBiasOnlyOnce(true) {
    fBremSplittingOperation = new GammaBiasingOperation("BremSplittingOperation");
}


void GammaBiasingOperator::StartRun() {
    fBremSplittingOperation->SetSplittingFactor(fSplittingFactor);
    G4cout << GetName() << " : starting run with brem. splitting factor = "
           << fSplittingFactor;
    if (fBiasPrimaryOnly) G4cout << ", biasing only primaries ";
    else
        G4cout << ", biasing primary and secondary tracks ";
    if (fBiasOnlyOnce) G4cout << ", biasing only once per track ";
    else
        G4cout << ", biasing several times per track ";
    G4cout << " . " << G4endl;
}

void GammaBiasingOperator::StartTracking(const G4Track * /* track */ ) {
    // -- reset the number of times the brem. splitting was applied:
    fNInteractions = 0;
}

G4VBiasingOperation *
GammaBiasingOperator::ProposeFinalStateBiasingOperation(const G4Track *track,
                                                        const G4BiasingProcessInterface * /* callingProcess */) {

    // -- Check if biasing of primary particle only is requested. If so, and
    // -- if particle is not a primary one, don't ask for biasing:
    if (fBiasPrimaryOnly && (track->GetParentID() != 0)) return nullptr;
    // -- Check if brem. splitting should be applied only once to the track,
    // -- and if so, and if brem. splitting already occured, don't ask for biasing:
    if (fBiasOnlyOnce && (fNInteractions > 0)) return nullptr;

    // -- Count the number of times the brem. splitting is applied:
    fNInteractions++;
    // -- Return the brem. splitting operation:
    return fBremSplittingOperation;
}