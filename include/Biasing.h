
#ifndef REST_BIASING_H
#define REST_BIASING_H

#include <TVector3.h>

#include <G4BiasingProcessInterface.hh>
#include <G4ParticleChange.hh>
#include <G4VBiasingOperation.hh>
#include <G4VBiasingOperator.hh>

class GammaBiasingOperation : public G4VBiasingOperation {
public:
    GammaBiasingOperation(const G4String &name, G4int splittingFactor, const TVector3 &biasingCenter);

    virtual ~GammaBiasingOperation();

public:
    virtual const G4VBiasingInteractionLaw *ProvideOccurenceBiasingInteractionLaw(
            const G4BiasingProcessInterface *, G4ForceCondition &) {
        return 0;
    }

    virtual G4VParticleChange *ApplyFinalStateBiasing(const G4BiasingProcessInterface *, const G4Track *,
                                                      const G4Step *, G4bool &);

    virtual G4double DistanceToApplyOperation(const G4Track *, G4double, G4ForceCondition *) { return DBL_MAX; }

    virtual G4VParticleChange *GenerateBiasingFinalState(const G4Track *, const G4Step *) { return 0; }

public:
    void SetSplittingFactor(G4int splittingFactor) { fSplittingFactor = splittingFactor; }

    G4int GetSplittingFactor() const { return fSplittingFactor; }

private:
    G4int fSplittingFactor;
    G4ParticleChange fParticleChange;
    TVector3 fBiasingCenter;
};


class GammaBiasingOperator : public G4VBiasingOperator {
public:
    GammaBiasingOperator(int splittingFactor, bool biasOnlyOnce, const TVector3 &biasingCenter);

    virtual ~GammaBiasingOperator() {}

public:
    virtual void StartRun();

    virtual void StartTracking(const G4Track *track);

    TVector3 GetBiasingCenter() const { return fBiasingCenter; }

private:
    virtual G4VBiasingOperation *ProposeNonPhysicsBiasingOperation(const G4Track *,
                                                                   const G4BiasingProcessInterface *) {
        return 0;
    }

    virtual G4VBiasingOperation *ProposeOccurenceBiasingOperation(const G4Track *,
                                                                  const G4BiasingProcessInterface *) {
        return 0;
    }

    virtual G4VBiasingOperation *ProposeFinalStateBiasingOperation(
            const G4Track *track, const G4BiasingProcessInterface *callingProcess);

private:
    using G4VBiasingOperator::OperationApplied;

private:
    GammaBiasingOperation *fBremSplittingOperation;
    G4int fSplittingFactor;
    G4bool fBiasOnlyOnce;
    G4int fNInteractions;
    TVector3 fBiasingCenter;
};

#endif //REST_BIASING_H
