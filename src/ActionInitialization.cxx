//
// Created by lobis on 6/14/2022.
//

#include "ActionInitialization.h"

#include "DetectorConstruction.h"
#include "EventAction.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "SimulationManager.h"
// #include "StackingAction.h"
#include "SteppingAction.h"
// #include "SteppingVerbose.h"
#include "G4RunManager.hh"
#include "TrackingAction.h"

using namespace std;

ActionInitialization::ActionInitialization(SimulationManager* simulationManager)
    : G4VUserActionInitialization(), fSimulationManager(simulationManager) {}

ActionInitialization::~ActionInitialization() = default;

void ActionInitialization::BuildForMaster() const { SetUserAction(new RunAction); }

void ActionInitialization::Build() const {
    auto detector = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    auto primaryGenerator = new PrimaryGeneratorAction(detector);

    auto restG4Metadata = fSimulationManager->fRestGeant4Metadata;

    if (restG4Metadata->GetParticleSource(0)->GetEnergyDistType() == "TH1D") {
        TString fileFullPath = (TString)restG4Metadata->GetParticleSource(0)->GetSpectrumFilename();

        TFile fin(fileFullPath);

        TString sptName = restG4Metadata->GetParticleSource(0)->GetSpectrumName();

        TH1D* h = (TH1D*)fin.Get(sptName);

        if (!h) {
            cout << "REST ERROR  when trying to find energy spectrum" << endl;
            cout << "File : " << fileFullPath << endl;
            cout << "Spectrum name : " << sptName << endl;
            exit(1);
        }

        fSimulationManager->initialEnergySpectrum = *h;

        Double_t minEnergy = restG4Metadata->GetParticleSource(0)->GetMinEnergy();
        if (minEnergy < 0) minEnergy = 0;

        Double_t maxEnergy = restG4Metadata->GetParticleSource(0)->GetMaxEnergy();
        if (maxEnergy < 0) maxEnergy = 0;

        // We set the initial spectrum energy provided from TH1D
        primaryGenerator->SetSpectrum(&(fSimulationManager->initialEnergySpectrum), minEnergy, maxEnergy);
    }

    if (restG4Metadata->GetParticleSource(0)->GetAngularDistType() == "TH1D") {
        TString fileFullPath = (TString)restG4Metadata->GetParticleSource(0)->GetAngularFilename();

        TFile fin(fileFullPath);

        TString sptName = restG4Metadata->GetParticleSource(0)->GetAngularName();
        TH1D* h = (TH1D*)fin.Get(sptName);

        if (!h) {
            cout << "REST ERROR  when trying to find angular spectrum" << endl;
            cout << "File : " << fileFullPath << endl;
            cout << "Spectrum name : " << sptName << endl;
            exit(1);
        }

        fSimulationManager->initialAngularDistribution = *h;

        // We set the initial angular distribution provided from TH1D
        primaryGenerator->SetAngularDistribution(&(fSimulationManager->initialAngularDistribution));
    }

    auto runAction = new RunAction();
    auto eventAction = new EventAction();
    auto trackingAction = new TrackingAction(runAction, eventAction);

    SetUserAction(new PrimaryGeneratorAction(detector));
    SetUserAction(runAction);
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction);
    // SetUserAction(new StackingAction);
    SetUserAction(trackingAction);

    /*
    G4EventManager::GetEventManager()->SetNumberOfAdditionalWaitingStacks(1);  // optical stack
     */
}

// G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const { return new SteppingVerbose; }