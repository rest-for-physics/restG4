
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

void Application::Run(const CommandLineParameters& commandLineParameters) {
    delete fSimulationManager;
    fSimulationManager = new SimulationManager();

    Bool_t saveAllEvents;

    Int_t nEvents;

    auto timeStart = chrono::steady_clock::now();

    const auto originalDirectory = filesystem::current_path();

    cout << "Current working directory: " << originalDirectory << endl;

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

    fSimulationManager->fRestGeant4Metadata = new TRestGeant4Metadata(inputRmlClean.c_str());
    fSimulationManager->fRestGeant4Metadata->SetGeant4Version(TRestTools::Execute("geant4-config --version"));

    if (!commandLineParameters.geometryFile.IsNull()) {
        fSimulationManager->fRestGeant4Metadata->SetGdmlFilename(commandLineParameters.geometryFile.Data());
    }

    // We need to process and generate a new GDML for several reasons.
    // 1. ROOT6 has problem loading math expressions in gdml file
    // 2. We allow file entities to be http remote files
    // 3. We retrieve the GDML and materials versions and associate to the
    // corresponding TRestGeant4Metadata members
    // 4. We support the use of system variables ${}
    auto gdml = new TRestGDMLParser();

    // This call will generate a new single file GDML output
    gdml->Load((string)fSimulationManager->fRestGeant4Metadata->GetGdmlFilename());

    // We redefine the value of the GDML file to be used in DetectorConstructor.
    fSimulationManager->fRestGeant4Metadata->SetGdmlFilename(gdml->GetOutputGDMLFile());
    fSimulationManager->fRestGeant4Metadata->SetGeometryPath("");

    fSimulationManager->fRestGeant4Metadata->SetGdmlReference(gdml->GetGDMLVersion());
    fSimulationManager->fRestGeant4Metadata->SetMaterialsReference(gdml->GetEntityVersion("materials"));

    fSimulationManager->fRestGeant4PhysicsLists = new TRestGeant4PhysicsLists(inputRmlClean.c_str());

    fSimulationManager->fRestRun = new TRestRun();
    fSimulationManager->fRestRun->LoadConfigFromFile(inputRmlClean);

    if (!commandLineParameters.outputFile.IsNull()) {
        fSimulationManager->fRestRun->SetOutputFileName(commandLineParameters.outputFile.Data());
    }

    filesystem::current_path(originalDirectory);

    TString runTag = fSimulationManager->fRestRun->GetRunTag();
    if (runTag == "Null" || runTag == "")
        fSimulationManager->fRestRun->SetRunTag(fSimulationManager->fRestGeant4Metadata->GetTitle());

    fSimulationManager->fRestRun->SetRunType("restG4");

    fSimulationManager->fRestRun->AddMetadata(fSimulationManager->fRestGeant4Metadata);
    fSimulationManager->fRestRun->AddMetadata(fSimulationManager->fRestGeant4PhysicsLists);
    fSimulationManager->fRestRun->PrintMetadata();

    fSimulationManager->fRestRun->FormOutputFile();

    fSimulationManager->fRestGeant4Event = new TRestGeant4Event();
    fSimulationManager->fRestGeant4SubEvent = new TRestGeant4Event();
    fSimulationManager->fRestRun->AddEventBranch(&fSimulationManager->fEvent);

    fSimulationManager->fRestGeant4Track = new TRestGeant4Track();

    // choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    long seed = fSimulationManager->fRestGeant4Metadata->GetSeed();
    CLHEP::HepRandom::setTheSeed(seed);

    auto runManager = new G4RunManager;

    auto detector = new DetectorConstruction(fSimulationManager);

    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList(fSimulationManager->fRestGeant4PhysicsLists));

    runManager->SetUserInitialization(new ActionInitialization(fSimulationManager));

    auto step = (SteppingAction*)G4RunManager::GetRunManager()->GetUserSteppingAction();

    auto primaryGenerator =
        (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif

    nEvents = fSimulationManager->fRestGeant4Metadata->GetNumberOfEvents();

    time_t systime = time(nullptr);
    fSimulationManager->fRestRun->SetStartTimeStamp((Double_t)systime);

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

        fSimulationManager->fRestRun->GetOutputFile()->cd();

    }

    else if (nEvents == 0)  // define visualization and UI terminal for interactive mode
    {
#ifdef G4UI_USE
        cout << "Entering vis mode.." << endl;
        auto ui = new G4UIExecutive(commandLineParameters.cmdArgc, commandLineParameters.cmdArgv);
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
    fSimulationManager->fRestRun->GetOutputFile()->cd();

#ifdef G4VIS_USE
    delete visManager;
#endif

    // job termination
    delete runManager;

    systime = time(nullptr);
    fSimulationManager->fRestRun->SetEndTimeStamp((Double_t)systime);
    const TString filename =
        TRestTools::ToAbsoluteName(fSimulationManager->fRestRun->GetOutputFileName().Data());
    fSimulationManager->fRestRun->UpdateOutputFile();
    fSimulationManager->fRestRun->CloseFile();
    fSimulationManager->fRestRun->PrintMetadata();

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
    auto timeEnd = chrono::steady_clock::now();
    cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(timeEnd - timeStart).count()
         << " seconds" << endl;
}
