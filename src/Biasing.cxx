
#include "Biasing.h"

#include <G4ParticleChangeForLoss.hh>

// Operation

GammaBiasingOperation::GammaBiasingOperation(const G4String &name, int splittingFactor,
                                             const TVector3 &biasingCenter)
        : G4VBiasingOperation(name),
          fSplittingFactor(splittingFactor),
          fParticleChange(),
          fBiasingCenter(biasingCenter) {}

GammaBiasingOperation::~GammaBiasingOperation() = default;

G4VParticleChange *GammaBiasingOperation::ApplyFinalStateBiasing(
        const G4BiasingProcessInterface *callingProcess, const G4Track *track, const G4Step *step, G4bool &) {
    G4VParticleChange *processFinalState = callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);

    if (fSplittingFactor == 1) return processFinalState;

    // special case: no secondaries
    if (processFinalState->GetNumberOfSecondaries() == 0) return processFinalState;

    auto actualParticleChange = (G4ParticleChangeForLoss *) processFinalState;

    fParticleChange.Initialize(*track);

    // -- Store electron final state:
    fParticleChange.ProposeTrackStatus(actualParticleChange->GetTrackStatus());
    fParticleChange.ProposeEnergy(actualParticleChange->GetProposedKineticEnergy());
    fParticleChange.ProposeMomentumDirection(actualParticleChange->GetProposedMomentumDirection());

    fParticleChange.SetSecondaryWeightByProcess(true);

    G4Track *gammaTrack = actualParticleChange->GetSecondary(0);

    int nSecondaries = 1;
    double weightFactor = 1.0;

    const auto diff =
            gammaTrack->GetPosition() - G4ThreeVector(fBiasingCenter.X(), fBiasingCenter.Y(), fBiasingCenter.Z());
    if (gammaTrack->GetMomentumDirection().dot(diff) < 0) {
        // pointing towards point of interest
        nSecondaries = fSplittingFactor;
        weightFactor = 1.0 / fSplittingFactor;
    } else {
        // pointing away from point of interest
        // random number between 0 and 1
        const double rand = G4UniformRand();
        // if random number is less than 1 / fSplittingFactor, keep alive
        weightFactor = fSplittingFactor; // increase weight
        if (rand < 1.0 / fSplittingFactor) {
            nSecondaries = 1;
        } else {
            nSecondaries = 0;
        }
    }

    G4double gammaWeight = track->GetWeight() * weightFactor;
    // G4cout << "Gamma weight: " << gammaWeight << " nSecondaries: " << nSecondaries << G4endl;
    gammaTrack->SetWeight(gammaWeight);

    fParticleChange.SetNumberOfSecondaries(nSecondaries);
    if (nSecondaries > 0) {
        fParticleChange.AddSecondary(gammaTrack);
    }

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

// Operator

#include "G4BiasingProcessInterface.hh"

GammaBiasingOperator::GammaBiasingOperator(int splittingFactor, bool biasOnlyOnce,
                                           const TVector3 &biasingCenter)
        : G4VBiasingOperator("BremSplittingOperator"),
          fSplittingFactor(splittingFactor),
          fBiasOnlyOnce(biasOnlyOnce),
          fBiasingCenter(biasingCenter) {
    fBremSplittingOperation =
            new GammaBiasingOperation("BremSplittingOperation", splittingFactor, biasingCenter);
}

void GammaBiasingOperator::StartRun() { fBremSplittingOperation->SetSplittingFactor(fSplittingFactor); }

void GammaBiasingOperator::StartTracking(const G4Track *) { fNInteractions = 0; }

G4VBiasingOperation *GammaBiasingOperator::ProposeFinalStateBiasingOperation(
        const G4Track *, const G4BiasingProcessInterface *) {
    if (fBiasOnlyOnce && (fNInteractions > 0)) return nullptr;

    fNInteractions++;

    return fBremSplittingOperation;
}
