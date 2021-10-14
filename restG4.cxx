
#include <TGeoVolume.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRestGDMLParser.h>
#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

#include <G4RunManager.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>

#include "ActionInitialization.h"
#include "CommandLineSetup.h"
#include "DetectorConstruction.h"
#include "EventAction.h"
#include "GlobalManager.h"
#include "PhysicsList.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "SteppingAction.h"
#include "SteppingVerbose.h"
#include "TrackingAction.h"
#include "spdlog/spdlog.h"

using namespace std;

// We define rest objects that will be used in Geant4
// TRestRun* restRun;
TRestGeant4Track* restTrack;
TRestGeant4Event *restG4Event, *subRestG4Event;
// TRestGeant4Metadata* restG4Metadata;
// TRestGeant4PhysicsLists* restPhysList;

Bool_t saveAllEvents;

const Int_t maxBiasingVolumes = 50;
Int_t biasing = 0;

// these histograms would be better placed inside TRestGeant4BiasingVolume
TH1D* biasingSpectrum[maxBiasingVolumes];
TH1D* angularDistribution[maxBiasingVolumes];
TH2D* spatialDistribution[maxBiasingVolumes];

TH1D initialEnergySpectrum;
TH1D initialAngularDistribution;

Int_t numberOfEvents;

int main(int argc, char** argv) {
    auto start_time = chrono::steady_clock::now();

    CommandLineParameters commandLineParameters = CommandLineSetup::ProcessParameters(argc, argv);
    CommandLineSetup::Print(commandLineParameters);

    GlobalManager::Instance()->InitializeFromConfigFile(commandLineParameters.rmlFile);

    auto restG4Metadata = GlobalManager::Instance()->GetRestGeant4Metadata();
    auto restRun = GlobalManager::Instance()->GetRestRun();
    // We need to process and generate a new GDML for several reasons.
    // 1. ROOT6 has problem loading math expressions in gdml file
    // 2. We allow file entities to be http remote files
    // 3. We retrieve the GDML and materials versions and associate to the
    // corresponding TRestGeant4Metadata members
    // 4. We support the use of system variables ${}
    auto gdml = new TRestGDMLParser();

    // This call will generate a new single file GDML output
    gdml->Load((string)GlobalManager::Instance()->GetRestGeant4Metadata()->Get_GDML_Filename());

    restTrack = new TRestGeant4Track();

    biasing = restG4Metadata->GetNumberOfBiasingVolumes();
    for (int i = 0; i < biasing; i++) {
        TString spctName = "Bias_Spectrum_" + TString(Form("%d", i));
        TString angDistName = "Bias_Angular_Distribution_" + TString(Form("%d", i));
        TString spatialDistName = "Bias_Spatial_Distribution_" + TString(Form("%d", i));

        Double_t maxEnergy = restG4Metadata->GetBiasingVolume(i).GetMaxEnergy();
        Double_t minEnergy = restG4Metadata->GetBiasingVolume(i).GetMinEnergy();
        Int_t nbins = (Int_t)(maxEnergy - minEnergy);

        Double_t biasSize = restG4Metadata->GetBiasingVolume(i).GetBiasingVolumeSize();
        TString biasType = restG4Metadata->GetBiasingVolume(i).GetBiasingVolumeType();

        cout << "Initializing biasing histogram : " << spctName << endl;
        biasingSpectrum[i] = new TH1D(spctName, "Biasing gamma spectrum", nbins, minEnergy, maxEnergy);
        angularDistribution[i] = new TH1D(angDistName, "Biasing angular distribution", 150, 0, M_PI / 2);

        if (biasType == "virtualSphere")
            spatialDistribution[i] =
                new TH2D(spatialDistName, "Biasing spatial (virtualSphere) distribution ", 100, -M_PI, M_PI,
                         100, 0, M_PI);
        else if (biasType == "virtualBox")
            spatialDistribution[i] =
                new TH2D(spatialDistName, "Biasing spatial (virtualBox) distribution", 100, -biasSize / 2.,
                         biasSize / 2., 100, -biasSize / 2., biasSize / 2.);
        else
            spatialDistribution[i] =
                new TH2D(spatialDistName, "Biasing spatial distribution", 100, -1, 1, 100, -1, 1);
    }

    G4VSteppingVerbose::SetInstance(new SteppingVerbose);

    auto runManager = new G4RunManager();

    runManager->SetUserInitialization(new DetectorConstruction);
    runManager->SetUserInitialization(
        new PhysicsList(GlobalManager::Instance()->GetRestGeant4PhysicsLists()));

    runManager->SetUserInitialization(new ActionInitialization);

    auto primaryGeneratorAction =
        (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

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

        initialEnergySpectrum = *h;

        Double_t minEnergy = restG4Metadata->GetParticleSource(0)->GetMinEnergy();
        if (minEnergy < 0) minEnergy = 0;

        Double_t maxEnergy = restG4Metadata->GetParticleSource(0)->GetMaxEnergy();
        if (maxEnergy < 0) maxEnergy = 0;

        // We set the initial spectrum energy provided from TH1D
        primaryGeneratorAction->SetSpectrum(&initialEnergySpectrum, minEnergy, maxEnergy);
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

        initialAngularDistribution = *h;

        // We set the initial angular distribution provided from TH1D
        primaryGeneratorAction->SetAngularDistribution(&initialAngularDistribution);
    }

    G4UImanager* UIManager = G4UImanager::GetUIpointer();
    G4VisManager* visManager = nullptr;
    G4UIExecutive* ui = nullptr;

    spdlog::info("Initializing run");
    UIManager->ApplyCommand("/run/initialize");

    if (commandLineParameters.interactive) {
        GlobalManager::Instance()->SetSaveFlag(false);
        GlobalManager::Instance()->SetInteractiveFlag(true);
        spdlog::info("Initializing G4VisExecutive");
        visManager = new G4VisExecutive();
        visManager->Initialize();

        ui = new G4UIExecutive(argc, argv);
    }
    if (visManager) {
        spdlog::info("Running visualization macros");
        CommandLineSetup::RunVisMacro();
        if (ui->IsGUI()) CommandLineSetup::RunGUIMacro();
        ui->SessionStart();
        delete ui;
    }

    // G4VSteppingVerbose::SetInstance(new SteppingVerbose);

    // runManager->Initialize();

    numberOfEvents = restG4Metadata->GetNumberOfEvents();
    // We pass the volume definition to Stepping action so that it records gammas
    // entering in We pass also the biasing spectrum so that gammas energies
    // entering the volume are recorded
    /*
    if (biasing) {
        step->SetBiasingVolume(restG4Metadata->GetBiasingVolume(biasing - 1));
        step->SetBiasingSpectrum(biasingSpectrum[biasing - 1]);
        step->SetAngularDistribution(angularDistribution[biasing - 1]);
        step->SetSpatialDistribution(spatialDistribution[biasing - 1]);
    }
*/
    time_t systime = time(nullptr);
    restRun->SetStartTimeStamp((Double_t)systime);

    if (!commandLineParameters.interactive) {
        spdlog::info("Starting batch mode with {} events to be processed", numberOfEvents);

        UIManager->ApplyCommand("/tracking/verbose 0");
        UIManager->ApplyCommand("/run/beamOn " + std::to_string(numberOfEvents));

        restRun->GetOutputFile()->cd();

        /*
        if (biasing) {
            cout << "Biasing id: " << biasing - 1 << endl;
            step->GetBiasingVolume().PrintBiasingVolume();
            cout << "Number of events that reached the biasing volume : "
                 << (Int_t)(biasingSpectrum[biasing - 1]->Integral()) << endl;
            cout << endl;
            cout << endl;
            biasing--;
        }
        while (biasing) {
            restG4Metadata->RemoveParticleSources();

            auto src = new TRestGeant4ParticleSource();
            src->SetParticleName("gamma");
            src->SetEnergyDistType("TH1D");
            src->SetAngularDistType("TH1D");
            restG4Metadata->AddParticleSource(src);

            // We set the spectrum from previous biasing volume inside the primary generator
            primaryGeneratorAction->SetSpectrum(biasingSpectrum[biasing]);
            // And we set the angular distribution
            primaryGeneratorAction->SetAngularDistribution(angularDistribution[biasing]);

            // We re-define the generator inside restG4Metadata to be launched from the biasing volume
            restG4Metadata->SetGeneratorType(
                restG4Metadata->GetBiasingVolume(biasing).GetBiasingVolumeType());
            double size = restG4Metadata->GetBiasingVolume(biasing).GetBiasingVolumeSize();
            restG4Metadata->SetGeneratorSize(TVector3(size, size, size));
            // restG4Metadata->GetBiasingVolume( biasing-1 ).PrintBiasingVolume();

            // Defining biasing the number of event to be re-launched
            Double_t biasingFactor = restG4Metadata->GetBiasingVolume(biasing - 1).GetBiasingFactor();
            cout << "Biasing id: " << biasing - 1 << ", Events to be launched : "
                 << (Int_t)(biasingSpectrum[biasing]->Integral() * biasingFactor) << endl;

            sprintf(tmp, "/run/beamOn %d", (Int_t)(biasingSpectrum[biasing]->Integral() * biasingFactor));
            command = tmp;

            restRun->GetOutputFile()->cd();

            biasingSpectrum[biasing]->Write();
            angularDistribution[biasing]->Write();
            spatialDistribution[biasing]->Write();

            // We pass the volume definition to Stepping action so that it records
            // gammas entering in We pass also the biasing spectrum so that gammas
            // energies entering the volume are recorded
            if (biasing) {
                step->SetBiasingVolume(restG4Metadata->GetBiasingVolume(biasing - 1));
                step->SetBiasingSpectrum(biasingSpectrum[biasing - 1]);
                step->SetAngularDistribution(angularDistribution[biasing - 1]);
                step->SetSpatialDistribution(spatialDistribution[biasing - 1]);
            }
            ui->ApplyCommand(command);
            step->GetBiasingVolume().PrintBiasingVolume();
            cout << "Number of events that reached the biasing volume : "
                 << (Int_t)(biasingSpectrum[biasing - 1]->Integral()) << endl;
            cout << endl;
            cout << endl;
            biasing--;
        }
        */
    }

    restRun->GetOutputFile()->cd();

    delete visManager;
    delete runManager;

    systime = time(nullptr);
    restRun->SetEndTimeStamp((Double_t)systime);
    TString Filename = restRun->GetOutputFileName();
    cout << "CLOSING: " << Filename << endl;
    restRun->UpdateOutputFile();
    restRun->CloseFile();
    restRun->PrintMetadata();
    delete restRun;

    ////////// Writing the geometry in TGeoManager format to the ROOT file
    ////////// Need to fork and do it in child process, to prevent occasional seg.fault
    pid_t pid;
    pid = fork();
    if (pid < 0) {
        perror("fork error:");
        exit(1);
    }
    // child process
    if (pid == 0) {
        // writing the geometry object
        freopen("/dev/null", "w", stdout);
        freopen("/dev/null", "w", stderr);
        Console::CompatibilityMode = true;

        // We wait the father process ends properly
        sleep(5);

        // Then we just add the geometry
        TFile* f1 = new TFile(Filename, "update");
        auto gdml = GlobalManager::Instance()->GetRestGDMLParser();
        TGeoManager* geo2 = gdml->CreateGeoM();

        f1->cd();
        geo2->SetName("Geometry");
        geo2->Write();
        f1->Close();
        exit(0);
    }
    // father process
    else {
        int stat_val = 0;
        pid_t child_pid;

        printf("Writing geometry ... \n");
    }

    cout << "============== Generated file: " << Filename << " ==============" << endl;
    auto end_time = chrono::steady_clock::now();
    cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(end_time - start_time).count()
         << " seconds" << endl;

    return 0;
}
