
#include "Application.h"

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
#include <G4UImanager.hh>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>

#include "ActionInitialization.h"
#include "CommandLineSetup.h"
#include "DetectorConstruction.h"
#include "EventAction.h"
#include "PhysicsList.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "SimulationManager.h"
#include "SteppingAction.h"
#include "TrackingAction.h"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

using namespace std;

void Application::Run(int argc, char** argv) {
    auto simulationManager = SimulationManager::Instance();

    TRestRun* restRun = simulationManager->fRestRun;
    TRestGeant4Track* restTrack = simulationManager->fRestGeant4Track;
    TRestGeant4Event* restG4Event = simulationManager->fRestGeant4Event;
    TRestGeant4Event* subRestG4Event = simulationManager->fRestGeant4SubEvent;
    TRestGeant4Metadata* restG4Metadata = simulationManager->fRestGeant4Metadata;
    TRestGeant4PhysicsLists* restPhysList = simulationManager->fRestGeant4PhysicsLists;
    Int_t& biasing = simulationManager->fBiasing;

    auto biasingSpectrum = simulationManager->biasingSpectrum;
    auto angularDistribution = simulationManager->biasingSpectrum;
    auto spatialDistribution = simulationManager->spatialDistribution;

    Bool_t saveAllEvents;

    Int_t nEvents;

    auto start_time = chrono::steady_clock::now();

    const auto originalDirectory = filesystem::current_path();

    cout << "Current working directory: " << originalDirectory << endl;

    CommandLineParameters commandLineParameters = CommandLineSetup::ProcessParameters(argc, argv);
    CommandLineSetup::Print(commandLineParameters);

    /// Separating relative path and pure RML filename
    char* inputConfigFile = const_cast<char*>(commandLineParameters.rmlFile.Data());

    if (!TRestTools::CheckFileIsAccessible(inputConfigFile)) {
        cout << "Input rml file: " << inputConfigFile << " not found, please check file name" << endl;
        exit(1);
    }

    const auto [inputRmlPath, inputRmlClean] = TRestTools::SeparatePathAndName(inputConfigFile);

    if (!filesystem::path(inputRmlPath).empty()) {
        filesystem::current_path(inputRmlPath);
    }

    restG4Metadata = new TRestGeant4Metadata(inputRmlClean.c_str());
    restG4Metadata->SetGeant4Version(TRestTools::Execute("geant4-config --version"));

    if (!commandLineParameters.geometryFile.IsNull()) {
        restG4Metadata->SetGdmlFilename(commandLineParameters.geometryFile.Data());
    }

    // We need to process and generate a new GDML for several reasons.
    // 1. ROOT6 has problem loading math expressions in gdml file
    // 2. We allow file entities to be http remote files
    // 3. We retrieve the GDML and materials versions and associate to the
    // corresponding TRestGeant4Metadata members
    // 4. We support the use of system variables ${}
    auto gdml = new TRestGDMLParser();

    // This call will generate a new single file GDML output
    gdml->Load((string)restG4Metadata->GetGdmlFilename());

    // We redefine the value of the GDML file to be used in DetectorConstructor.
    restG4Metadata->SetGdmlFilename(gdml->GetOutputGDMLFile());
    restG4Metadata->SetGeometryPath("");

    restG4Metadata->SetGdmlReference(gdml->GetGDMLVersion());
    restG4Metadata->SetMaterialsReference(gdml->GetEntityVersion("materials"));

    restPhysList = new TRestGeant4PhysicsLists(inputRmlClean.c_str());

    restRun = new TRestRun();
    restRun->LoadConfigFromFile(inputRmlClean);

    if (!commandLineParameters.outputFile.IsNull()) {
        restRun->SetOutputFileName(commandLineParameters.outputFile.Data());
    }

    filesystem::current_path(originalDirectory);

    TString runTag = restRun->GetRunTag();
    if (runTag == "Null" || runTag == "") restRun->SetRunTag(restG4Metadata->GetTitle());

    restRun->SetRunType("restG4");

    restRun->AddMetadata(restG4Metadata);
    restRun->AddMetadata(restPhysList);
    restRun->PrintMetadata();

    restRun->FormOutputFile();

    restG4Event = new TRestGeant4Event();
    subRestG4Event = new TRestGeant4Event();
    restRun->AddEventBranch(subRestG4Event);

    restTrack = new TRestGeant4Track();

    biasing = restG4Metadata->GetNumberOfBiasingVolumes();
    for (int i = 0; i < biasing; i++) {
        TString spectrumName = "Bias_Spectrum_" + TString(Form("%d", i));
        TString angDistName = "Bias_Angular_Distribution_" + TString(Form("%d", i));
        TString spatialDistName = "Bias_Spatial_Distribution_" + TString(Form("%d", i));

        Double_t maxEnergy = restG4Metadata->GetBiasingVolume(i).GetMaxEnergy();
        Double_t minEnergy = restG4Metadata->GetBiasingVolume(i).GetMinEnergy();
        auto nBins = (Int_t)(maxEnergy - minEnergy);

        Double_t biasSize = restG4Metadata->GetBiasingVolume(i).GetBiasingVolumeSize();
        TString biasType = restG4Metadata->GetBiasingVolume(i).GetBiasingVolumeType();

        cout << "Initializing biasing histogram : " << spectrumName << endl;
        biasingSpectrum[i] = new TH1D(spectrumName, "Biasing gamma spectrum", nBins, minEnergy, maxEnergy);
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

    // choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    long seed = restG4Metadata->GetSeed();
    CLHEP::HepRandom::setTheSeed(seed);

    auto runManager = new G4RunManager;

    auto detector = new DetectorConstruction();

    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList(restPhysList));

    runManager->SetUserInitialization(new ActionInitialization(SimulationManager::Instance()));

    auto step = (SteppingAction*)G4RunManager::GetRunManager()->GetUserSteppingAction();
    auto primaryGenerator =
        (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif

    nEvents = restG4Metadata->GetNumberOfEvents();
    // We pass the volume definition to Stepping action so that it records gammas
    // entering in We pass also the biasing spectrum so that gammas energies
    // entering the volume are recorded
    if (biasing) {
        step->SetBiasingVolume(restG4Metadata->GetBiasingVolume(biasing - 1));
        step->SetBiasingSpectrum(biasingSpectrum[biasing - 1]);
        step->SetAngularDistribution(angularDistribution[biasing - 1]);
        step->SetSpatialDistribution(spatialDistribution[biasing - 1]);
    }

    time_t systime = time(nullptr);
    restRun->SetStartTimeStamp((Double_t)systime);

    cout << "Number of events : " << nEvents << endl;
    if (nEvents > 0)  // batch mode
    {
        G4String command = "/tracking/verbose 0";
        UI->ApplyCommand(command);
        command = "/run/initialize";
        UI->ApplyCommand(command);

        char tmp[256];
        sprintf(tmp, "/run/beamOn %d", nEvents);

        command = tmp;
        UI->ApplyCommand(command);

        restRun->GetOutputFile()->cd();

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
            primaryGenerator->SetSpectrum(biasingSpectrum[biasing]);
            // And we set the angular distribution
            primaryGenerator->SetAngularDistribution(angularDistribution[biasing]);

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
            UI->ApplyCommand(command);
            step->GetBiasingVolume().PrintBiasingVolume();
            cout << "Number of events that reached the biasing volume : "
                 << (Int_t)(biasingSpectrum[biasing - 1]->Integral()) << endl;
            cout << endl;
            cout << endl;
            biasing--;
        }
    }

    else if (nEvents == 0)  // define visualization and UI terminal for interactive mode
    {
        cout << "Entering vis mode.." << endl;
#ifdef G4UI_USE
        auto ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        cout << "Executing G4 macro : /control/execute macros/vis.mac" << endl;
        UI->ApplyCommand("/control/execute macros/vis.mac");
#endif
        ui->SessionStart();
        delete ui;
#endif
    }

    else  // nEvents == -1
    {
        cout << "++++++++++ ERROR +++++++++" << endl;
        cout << "++++++++++ ERROR +++++++++" << endl;
        cout << "++++++++++ ERROR +++++++++" << endl;
        cout << "The number of events to be simulated was not recognized properly!" << endl;
        cout << "Make sure you did not forget the number of events entry in TRestGeant4Metadata." << endl;
        cout << endl;
        cout << " ... or the parameter is properly constructed/interpreted." << endl;
        cout << endl;
        cout << "It should be something like : " << endl;
        cout << endl;
        cout << R"( <parameter name ="nEvents" value="100"/>)" << endl;
        cout << "++++++++++ ERROR +++++++++" << endl;
        cout << "++++++++++ ERROR +++++++++" << endl;
        cout << "++++++++++ ERROR +++++++++" << endl;
        cout << endl;
    }
    restRun->GetOutputFile()->cd();

#ifdef G4VIS_USE
    delete visManager;
#endif

    // job termination
    delete runManager;

    systime = time(nullptr);
    restRun->SetEndTimeStamp((Double_t)systime);
    TString filename = TRestTools::ToAbsoluteName(restRun->GetOutputFileName().Data());

    restRun->UpdateOutputFile();
    restRun->CloseFile();
    restRun->PrintMetadata();
    delete restRun;

    delete restG4Event;
    delete subRestG4Event;
    delete restTrack;

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

        REST_Display_CompatibilityMode = true;

        // We wait the father process ends properly
        sleep(5);

        // Then we just add the geometry
        auto file = new TFile(filename, "update");
        TGeoManager* geoManager = gdml->CreateGeoManager();

        file->cd();
        geoManager->SetName("Geometry");
        geoManager->Write();
        file->Close();
        exit(0);
    }
    // father process
    else {
        int stat_val = 0;
        pid_t child_pid;

        printf("Writing geometry ... \n");
    }

    cout << "============== Generated file: " << filename << " ==============" << endl;
    auto end_time = chrono::steady_clock::now();
    cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(end_time - start_time).count()
         << " seconds" << endl;
}