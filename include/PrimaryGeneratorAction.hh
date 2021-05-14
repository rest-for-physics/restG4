#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <TH1D.h>

#include <fstream>
#include <iostream>

#include "DetectorConstruction.hh"
#include "G4IonTable.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

using namespace std;

const int nSpct = 3000;

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
   public:
    PrimaryGeneratorAction(DetectorConstruction* pDetector);
    ~PrimaryGeneratorAction();

   public:
    virtual void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() { return fParticleGun; };

    void SetSpectrum(TH1D* spt, double eMin = 0, double eMax = 0) {
        TString xLabel = (TString)spt->GetXaxis()->GetTitle();

        if (xLabel.Contains("MeV")) {
            energyFactor = 1.e3;
        } else if (xLabel.Contains("GeV")) {
            energyFactor = 1.e6;
        } else {
            energyFactor = 1.;
        }

        fSpectrum = spt;
        fSpectrumIntegral = fSpectrum->Integral();

        startEnergyBin = 1;
        endEnergyBin = fSpectrum->GetNbinsX();

        if (eMin > 0) {
            for (int i = startEnergyBin; i <= endEnergyBin; i++) {
                if (fSpectrum->GetBinCenter(i) > eMin) {
                    startEnergyBin = i;
                    break;
                }
            }
        }

        if (eMax > 0) {
            for (int i = startEnergyBin; i <= endEnergyBin; i++) {
                if (fSpectrum->GetBinCenter(i) > eMax) {
                    endEnergyBin = i;
                    break;
                }
            }
        }

        fSpectrumIntegral = fSpectrum->Integral(startEnergyBin, endEnergyBin);
    }

    void SetAngularDistribution(TH1D* ang) { fAngularDistribution = ang; }

   private:
    G4ParticleGun* fParticleGun;
    DetectorConstruction* fDetector;

    G4ParticleDefinition* fParticle = nullptr;

    TH1D* fSpectrum;
    TH1D* fAngularDistribution;

    Int_t startEnergyBin;
    Int_t endEnergyBin;
    Double_t fSpectrumIntegral;

    Int_t nBiasingVolumes;

    Double_t energyFactor;

    Double_t lastEnergy;

    void SetParticlePosition();
    void SetParticlePosition(int n);
    G4ParticleDefinition* SetParticleDefinition(int n);
    void SetParticleEnergy(int n);
    void SetParticleDirection(int n);

    G4ThreeVector GetIsotropicVector();
    Double_t GetAngle(G4ThreeVector x, G4ThreeVector y);
    Double_t GetCosineLowRandomThetaAngle();

    G4String fParType;
    G4String fGenType;
    G4double fParEnergy;
    G4double fParGenerator;

    G4String fSpctFilename;

    G4int gammaSpectrum[nSpct];

    G4int nCollections;
};

#endif
