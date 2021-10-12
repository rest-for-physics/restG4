
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <TF3.h>
#include <TH1D.h>

#include <G4IonTable.hh>
#include <G4ParticleGun.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <fstream>
#include <globals.hh>
#include <iostream>

#include "DetectorConstruction.h"
#include "TRestGeant4Particle.h"

using namespace std;

const int nSpct = 3000;

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
   public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

   public:
    virtual void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() { return fParticleGun; };

    void SetSpectrum(TH1D* spt, double eMin = 0, double eMax = 0);
    void SetGeneratorSpatialDensity(TString str);

    void SetAngularDistribution(TH1D* ang) { fAngularDistribution = ang; }

   private:
    TRestGeant4Metadata* fRestGeant4Metadata;

    vector<TRestGeant4Particle> fTempParticles;

    G4ParticleGun* fParticleGun;
    DetectorConstruction* fDetector;
    G4ParticleDefinition* fParticle = nullptr;

    TH1D* fSpectrum;
    TH1D* fAngularDistribution;
    TF3* fGeneratorSpatialDensityFunction;

    Int_t startEnergyBin;
    Int_t endEnergyBin;
    Double_t fSpectrumIntegral;

    Int_t nBiasingVolumes;

    Double_t energyFactor;

    Double_t lastEnergy;

    void SetParticlePosition();
    G4ParticleDefinition* SetParticleDefinition(Int_t particlesourceindex, TRestGeant4Particle p);
    void SetParticleEnergy(Int_t particlesourceindex, TRestGeant4Particle p);
    void SetParticleDirection(Int_t particlesourceindex, TRestGeant4Particle p);

    G4ThreeVector GetIsotropicVector();
    Double_t GetAngle(G4ThreeVector x, G4ThreeVector y);
    Double_t GetCosineLowRandomThetaAngle();

    void GenPositionOnGDMLVolume(double& x, double& y, double& z);
    void GenPositionOnGDMLSurface(double& x, double& y, double& z);
    void GenPositionOnBoxVolume(double& x, double& y, double& z);
    void GenPositionOnBoxSurface(double& x, double& y, double& z);
    void GenPositionOnSphereVolume(double& x, double& y, double& z);
    void GenPositionOnSphereSurface(double& x, double& y, double& z);
    void GenPositionOnCylinderVolume(double& x, double& y, double& z);
    void GenPositionOnCylinderSurface(double& x, double& y, double& z);
    void GenPositionOnPoint(double& x, double& y, double& z);
    void GenPositionOnWall(double& x, double& y, double& z);
    void GenPositionOnPlate(double& x, double& y, double& z);

    G4String fParType;
    G4String fGenType;
    G4double fParEnergy;
    G4double fParGenerator;

    G4String fSpctFilename;

    G4int gammaSpectrum[nSpct];

    G4int nCollections;
};

#endif
