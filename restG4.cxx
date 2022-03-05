// REST G4 main program. restG4.cxx
//
// Author : J. Galan
// Date : Jul-2015
//

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

#include "CommandLineSetup.h"
#include "DetectorConstruction.h"
#include "EventAction.h"
#include "PhysicsList.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "SteppingAction.h"
#include "TrackingAction.h"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

using namespace std;

// We define rest objects that will be used in Geant4
TRestRun* restRun;
TRestGeant4Track* restTrack;
TRestGeant4Event *restG4Event, *subRestG4Event;
TRestGeant4Metadata* restG4Metadata;
TRestGeant4PhysicsLists* restPhysList;

Bool_t saveAllEvents;

const Int_t maxBiasingVolumes = 50;
Int_t biasing = 0;

// This histograms would be better placed inside TRestGeant4BiasingVolume
TH1D* biasingSpectrum[maxBiasingVolumes];
TH1D* angularDistribution[maxBiasingVolumes];
TH2D* spatialDistribution[maxBiasingVolumes];

TH1D initialEnergySpectrum;
TH1D initialAngularDistribution;

Int_t N_events;

int main(int argc, char** argv) {
    auto start_time = chrono::steady_clock::now();

    char cwd[kMAXPATHLEN];
    cout << "Current working directory: " << getcwd(cwd, sizeof(cwd)) << endl;

    CommandLineParameters commandLineParameters = CommandLineSetup::ProcessParameters(argc, argv);
    CommandLineSetup::Print(commandLineParameters);

    /// Separating relative path and pure RML filename
    char* inputConfigFile = const_cast<char*>(commandLineParameters.rmlFile.Data());
    std::pair<string, string> pathAndRml = TRestTools::SeparatePathAndName(inputConfigFile);
    char* inputRMLClean = (char*)pathAndRml.second.data();

    TRestTools::ChangeDirectory(pathAndRml.first);

    restG4Metadata = new TRestGeant4Metadata(inputRMLClean);

    if (!commandLineParameters.geometryFile.IsNull()) {
        restG4Metadata->Set_GDML_Filename(commandLineParameters.geometryFile.Data());
    }

    string geant4Version = TRestTools::Execute("geant4-config --version");
    restG4Metadata->SetGeant4Version(geant4Version);

    // We need to process and generate a new GDML for several reasons.
    // 1. ROOT6 has problem loading math expressions in gdml file
    // 2. We allow file entities to be http remote files
    // 3. We retrieve the GDML and materials versions and associate to the
    // corresponding TRestGeant4Metadata members
    // 4. We support the use of system variables ${}
    auto gdml = new TRestGDMLParser();

    // This call will generate a new single file GDML output
    auto geometryFilename = restG4Metadata->Get_GDML_Filename();
    cout << "Geometry Filename: " << geometryFilename << endl;
    gdml->Load((string)restG4Metadata->Get_GDML_Filename());

    // We redefine the value of the GDML file to be used in DetectorConstructor.
    restG4Metadata->Set_GDML_Filename(gdml->GetOutputGDMLFile());
    restG4Metadata->SetGeometryPath("");

    restG4Metadata->Set_GDML_Reference(gdml->GetGDMLVersion());
    restG4Metadata->SetMaterialsReference(gdml->GetEntityVersion("materials"));

    restPhysList = new TRestGeant4PhysicsLists(inputRMLClean);

    restRun = new TRestRun();
    restRun->LoadConfigFromFile(inputRMLClean);

    if (!commandLineParameters.outputFile.IsNull()) {
        restRun->SetOutputFileName(commandLineParameters.outputFile.Data());
    }

    TRestTools::ReturnToPreviousDirectory();

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
    // }}}

    // choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    long seed = restG4Metadata->GetSeed();
    CLHEP::HepRandom::setTheSeed(seed);

    auto runManager = new G4RunManager;

    auto det = new DetectorConstruction();

    runManager->SetUserInitialization(det);
    runManager->SetUserInitialization(new PhysicsList(restPhysList));

    auto prim = new PrimaryGeneratorAction(det);

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
        prim->SetSpectrum(&initialEnergySpectrum, minEnergy, maxEnergy);
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
        prim->SetAngularDistribution(&initialAngularDistribution);
    }

    auto run = new RunAction(prim);

    auto event = new EventAction();

    auto track = new TrackingAction(run, event);
    auto step = new SteppingAction();

    runManager->SetUserAction(run);
    runManager->SetUserAction(prim);
    runManager->SetUserAction(event);
    runManager->SetUserAction(track);
    runManager->SetUserAction(step);

    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif

    N_events = restG4Metadata->GetNumberOfEvents();
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

    cout << "Events : " << N_events << endl;
    if (N_events > 0)  // batch mode
    {
        G4String command = "/tracking/verbose 0";
        UI->ApplyCommand(command);
        command = "/run/initialize";
        UI->ApplyCommand(command);

        char tmp[256];
        sprintf(tmp, "/run/beamOn %d", N_events);

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
            prim->SetSpectrum(biasingSpectrum[biasing]);
            // And we set the angular distribution
            prim->SetAngularDistribution(angularDistribution[biasing]);

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

    else if (N_events == 0)  // define visualization and UI terminal for interactive mode
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

    else  // N_events == -1
    {
        cout << "++++++++++ ERRORRRR +++++++++" << endl;
        cout << "++++++++++ ERRORRRR +++++++++" << endl;
        cout << "++++++++++ ERRORRRR +++++++++" << endl;
        cout << "The number of events to be simulated was not recongnized properly!" << endl;
        cout << "Make sure you did not forget the number of events entry in TRestGeant4Metadata." << endl;
        cout << endl;
        cout << " ... or the parameter is properly constructed/interpreted." << endl;
        cout << endl;
        cout << "It should be something like : " << endl;
        cout << endl;
        cout << " <parameter name =\"Nevents\" value=\"100\"/>" << endl;
        cout << "++++++++++ ERRORRRR +++++++++" << endl;
        cout << "++++++++++ ERRORRRR +++++++++" << endl;
        cout << "++++++++++ ERRORRRR +++++++++" << endl;
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
    TString Filename = restRun->GetOutputFileName();

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
        Console::CompatibilityMode = true;

        // We wait the father process ends properly
        sleep(5);

        // Then we just add the geometry
        TFile* f1 = new TFile(Filename, "update");
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
