
#ifndef REST_SIMULATIONMANAGER_H
#define REST_SIMULATIONMANAGER_H

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

class SimulationManager {
   public:
    SimulationManager();
    ~SimulationManager();

    TRestRun* fRestRun = nullptr;
    TRestGeant4Track* fRestGeant4Track = nullptr;
    TRestGeant4Event *fRestGeant4Event = nullptr, *fRestGeant4SubEvent = nullptr;
    TRestGeant4Metadata* fRestGeant4Metadata = nullptr;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists = nullptr;
    Int_t fBiasing = 0;

    static const Int_t maxBiasingVolumes = 50;

    // These histograms would be better placed inside TRestGeant4BiasingVolume
    TH1D* biasingSpectrum[maxBiasingVolumes];
    TH1D* angularDistribution[maxBiasingVolumes];
    TH2D* spatialDistribution[maxBiasingVolumes];

    TH1D initialEnergySpectrum;
    TH1D initialAngularDistribution;
};

#endif  // REST_SIMULATIONMANAGER_H
