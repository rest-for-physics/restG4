
#include <utility>

#include "GammaBiasingOperation.h"
#include "GammaBiasingOperator.h"

#include "G4ParticleChangeForLoss.hh"

GammaBiasingOperation::GammaBiasingOperation(G4String name)
        : G4VBiasingOperation(std::move(name)),
          fSplittingFactor(1),
          fParticleChange() {
}


GammaBiasingOperation::~GammaBiasingOperation() = default;

G4VParticleChange *GammaBiasingOperation::ApplyFinalStateBiasing(const G4BiasingProcessInterface *callingProcess,
                                                                 const G4Track *track,
                                                                 const G4Step *step,
                                                                 G4bool &) {

    G4VParticleChange *processFinalState =
            callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);

    if (fSplittingFactor == 1) return processFinalState;

    // special case: no secondaries
    if (processFinalState->GetNumberOfSecondaries() == 0) return processFinalState;

    auto actualParticleChange = (G4ParticleChangeForLoss *) processFinalState;

    fParticleChange.Initialize(*track);

    // -- Store electron final state:
    fParticleChange.
            ProposeTrackStatus(actualParticleChange->GetTrackStatus());
    fParticleChange.
            ProposeEnergy(actualParticleChange->GetProposedKineticEnergy());
    fParticleChange.
            ProposeMomentumDirection(actualParticleChange->GetProposedMomentumDirection());

    fParticleChange.SetSecondaryWeightByProcess(true);

    G4Track *gammaTrack = actualParticleChange->GetSecondary(0);


    // print gamma info (energy, direction)
    /*
    G4cout << "Gamma info. Weight: " << gammaTrack->GetWeight()
           << " Energy (keV): " << gammaTrack->GetKineticEnergy() / CLHEP::keV << " Direction: "
           << gammaTrack->GetMomentumDirection() << G4endl;
    */
    // if direction points towards (0,0,0)
    bool split = false;
    const auto diff = gammaTrack->GetPosition() - G4ThreeVector(0, 0, 0);
    if (gammaTrack->GetMomentumDirection().dot(diff) < 0) {
        // G4cout << "Gamma points towards (0,0,0)" << G4endl;
        split = true;
    }

    const int nSecondaries = split ? fSplittingFactor : 1;

    G4double gammaWeight = track->GetWeight() / nSecondaries;
    // G4cout << "Gamma weight: " << gammaWeight << " nSecondaries: " << nSecondaries << G4endl;
    gammaTrack->SetWeight(gammaWeight);


    fParticleChange.SetNumberOfSecondaries(nSecondaries);
    fParticleChange.AddSecondary(gammaTrack);

    actualParticleChange->Clear();

    G4int nCalls = 1;
    while (nCalls < nSecondaries) {
        processFinalState = callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);
        if (processFinalState->GetNumberOfSecondaries() == 1) {
            gammaTrack = processFinalState->GetSecondary(0);
            gammaTrack->SetWeight(gammaWeight);
            fParticleChange.AddSecondary(gammaTrack);
            nCalls++;
        }
            // -- very rare special case: we ignore for now.
        else if (processFinalState->GetNumberOfSecondaries() > 1) {
            for (G4int i = 0; i < processFinalState->GetNumberOfSecondaries(); i++)
                delete processFinalState->GetSecondary(i);
        }
        processFinalState->Clear();
    }

    return &fParticleChange;
}
