
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

using namespace std;

class TRestGeant4Metadata;
class OutputManager;

class SteppingAction : public G4UserSteppingAction {
   public:
    SteppingAction();
    ~SteppingAction();

    void SetBiasingVolume(TRestGeant4BiasingVolume biasVol) { restBiasingVolume = biasVol; }
    void SetBiasingSpectrum(TH1D* bSpectrum) { fBiasingSpectrum = bSpectrum; }
    void SetAngularDistribution(TH1D* aDist) { fAngularDistribution = aDist; }
    void SetSpatialDistribution(TH2D* sDist) { fSpatialDistribution = sDist; }
    void UserSteppingAction(const G4Step*);

    TRestGeant4BiasingVolume GetBiasingVolume() { return restBiasingVolume; }
    TH1D* GetBiasingSpectrum() { return fBiasingSpectrum; }
    TH1D* GetAngularDistribution() { return fAngularDistribution; }
    TH2D* GetSpatialDistribution() { return fSpatialDistribution; }

    G4TrackVector* GetfSecondary();

   private:
    TRestGeant4Metadata* fRestGeant4Metadata;
    OutputManager* fOutputManager;

    G4String nom_vol, nom_part, nom_proc;
    G4double dif_ener, ener_dep, ener, eKin;
    G4int trackID, parentID;

    G4ThreeVector previousDirection;
    TH1D* fBiasingSpectrum;
    TH1D* fAngularDistribution;
    TH2D* fSpatialDistribution;

    Double_t absDouble(Double_t x) {
        if (x < 0) return -x;
        return x;
    }

    TRestGeant4BiasingVolume restBiasingVolume;
};
#endif
