
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <TF3.h>
#include <TH1D.h>
#include <TRandom.h>
#include <TRestGeant4Particle.h>

#include <G4IonTable.hh>
#include <G4ParticleGun.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <fstream>
#include <globals.hh>
#include <iostream>
#include <mutex>

#include "DetectorConstruction.h"

class G4Event;
class SimulationManager;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
   public:
    PrimaryGeneratorAction(SimulationManager*);
    ~PrimaryGeneratorAction();

   public:
    virtual void GeneratePrimaries(G4Event*);

    void SetEnergyDistributionHistogram(const TH1D* h, double eMin = 0, double eMax = 0);
    inline void SetAngularDistributionHistogram(const TH1D* h) { fAngularDistributionHistogram = h; }

    void SetGeneratorSpatialDensity(TString str);

   private:
    SimulationManager* fSimulationManager;
    std::mutex fMutex;

    std::vector<TRestGeant4Particle> fTempParticles;

    G4ParticleGun fParticleGun;
    G4ParticleDefinition* fParticle = nullptr;

    const TH1D* fEnergyDistributionHistogram = nullptr;
    const TH1D* fAngularDistributionHistogram = nullptr;

    TF1* fEnergyDistributionFunction = nullptr;
    TF1* fAngularDistributionFunction = nullptr;
    TF2* fEnergyAndAngularDistributionFunction = nullptr;

    TF3* fGeneratorSpatialDensityFunction;

    Int_t startEnergyBin;
    Int_t endEnergyBin;
    Double_t fSpectrumIntegral;

    Double_t energyFactor;

    Double_t lastEnergy;

    TRandom* fRandom = nullptr;

    void SetParticlePosition();
    G4ParticleDefinition* SetParticleDefinition(Int_t particleSourceIndex,
                                                const TRestGeant4Particle& particle);
    void SetParticleEnergy(Int_t particleSourceIndex, const TRestGeant4Particle& particle);
    void SetParticleDirection(Int_t particleSourceIndex, const TRestGeant4Particle& particle);

    void SetParticleEnergyAndDirection(Int_t particleSourceIndex, const TRestGeant4Particle& particle);

    G4ThreeVector GetIsotropicVector() const;
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
    void GenPositionOnDisk(double& x, double& y, double& z);

    G4String fParType;
    G4String fGenType;
};

#endif
