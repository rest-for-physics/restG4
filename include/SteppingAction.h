
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include <TH1D.h>
#include <TH2D.h>
#include <TRestGeant4BiasingVolume.h>

#include <G4RunManager.hh>
#include <G4UserSteppingAction.hh>
#include <fstream>
#include <globals.hh>
#include <iostream>

class SteppingAction : public G4UserSteppingAction {
   public:
    SteppingAction();
    ~SteppingAction();

    void SetBiasingVolume(const TRestGeant4BiasingVolume& biasVol) { restBiasingVolume = biasVol; }
    void SetBiasingSpectrum(TH1D* bSpectrum) { fBiasingSpectrum = bSpectrum; }
    void SetAngularDistribution(TH1D* aDist) { fAngularDistribution = aDist; }
    void SetSpatialDistribution(TH2D* sDist) { fSpatialDistribution = sDist; }
    void UserSteppingAction(const G4Step*) override;

    TRestGeant4BiasingVolume GetBiasingVolume() { return restBiasingVolume; }
    TH1D* GetBiasingSpectrum() { return fBiasingSpectrum; }
    TH1D* GetAngularDistribution() { return fAngularDistribution; }
    TH2D* GetSpatialDistribution() { return fSpatialDistribution; }

   private:
    G4double dif_ener, ener_dep, ener, eKin;
    G4int trackID, parentID;

    TH1D* fBiasingSpectrum;
    TH1D* fAngularDistribution;
    TH2D* fSpatialDistribution;

    static Double_t absDouble(Double_t x) {
        if (x < 0) return -x;
        return x;
    }

    TRestGeant4BiasingVolume restBiasingVolume;

    G4ThreeVector previousDirection;
};
#endif
