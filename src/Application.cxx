
#include "Application.h"

#include <TGeoVolume.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TRestGDMLParser.h>
#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

#include <G4RunManager.hh>
#ifndef GEANT4_WITHOUT_G4RunManagerFactory
#include <G4RunManagerFactory.hh>
#endif
#include <G4UImanager.hh>
#include <G4VSteppingVerbose.hh>
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
#include "SteppingVerbose.h"
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

    const auto timeStart = chrono::steady_clock::now();

    const auto originalDirectory = filesystem::current_path();

    cout << "Current working directory: " << originalDirectory << endl;

    CommandLineSetup::Print(commandLineParameters);

    /// Separating relative path and pure RML filename
    const char* inputConfigFile = const_cast<char*>(commandLineParameters.rmlFile.Data());

    if (!TRestTools::CheckFileIsAccessible(inputConfigFile)) {
        cout << "Input rml file: " << inputConfigFile << " not found, please check file name" << endl;
        exit(1);
    }

    const auto [inputRmlPath, inputRmlClean] = TRestTools::SeparatePathAndName(inputConfigFile);

    if (!filesystem::path(inputRmlPath).empty()) {
        filesystem::current_path(inputRmlPath);
    }

    auto metadata = new TRestGeant4Metadata(inputRmlClean.c_str());
    fSimulationManager->SetRestMetadata(metadata);

    metadata->SetGeant4Version(TRestTools::Execute("geant4-config --version"));

    if (!commandLineParameters.geometryFile.IsNull()) {
        metadata->SetGdmlFilename(commandLineParameters.geometryFile.Data());
    }

    // We need to process and generate a new GDML for several reasons.
    // 1. ROOT6 has problem loading math expressions in gdml file
    // 2. We allow file entities to be http remote files
    // 3. We retrieve the GDML and materials versions and associate to the
    // corresponding TRestGeant4Metadata members
    // 4. We support the use of system variables ${}
    auto gdml = new TRestGDMLParser();

    // This call will generate a new single file GDML output
    gdml->Load((string)metadata->GetGdmlFilename());

    // We redefine the value of the GDML file to be used in DetectorConstructor.
    metadata->SetGdmlFilename(gdml->GetOutputGDMLFile());
    metadata->SetGeometryPath("");

    metadata->SetGdmlReference(gdml->GetGDMLVersion());
    metadata->SetMaterialsReference(gdml->GetEntityVersion("materials"));

    auto physicsLists = new TRestGeant4PhysicsLists(inputRmlClean.c_str());
    fSimulationManager->SetRestPhysicsLists(physicsLists);

    auto run = new TRestRun();
    fSimulationManager->SetRestRun(run);

    run->LoadConfigFromFile(inputRmlClean);

    if (!commandLineParameters.outputFile.IsNull()) {
        run->SetOutputFileName(commandLineParameters.outputFile.Data());
    }

    filesystem::current_path(originalDirectory);

    TString runTag = run->GetRunTag();
    if (runTag == "Null" || runTag == "") {
        run->SetRunTag(metadata->GetTitle());
    }

    run->SetRunType("restG4");

    run->AddMetadata(fSimulationManager->GetRestMetadata());
    run->AddMetadata(fSimulationManager->GetRestPhysicsLists());
    run->PrintMetadata();

    run->FormOutputFile();

    run->GetOutputFile()->cd();

    cout << "Writing geometry into output file" << endl;
    TGeoManager* geoManager = gdml->CreateGeoManager();
    if (geoManager) {
        geoManager->Write("Geometry");
        // delete geoManager;
    } else {
        cout << "Error: could not save geometry" << endl;
        exit(1);
    }

    run->AddEventBranch(&fSimulationManager->fEvent);

    // choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    long seed = metadata->GetSeed();
    CLHEP::HepRandom::setTheSeed(seed);

    G4VSteppingVerbose::SetInstance(new SteppingVerbose(fSimulationManager));

#ifndef GEANT4_WITHOUT_G4RunManagerFactory
    auto runManagerType = G4RunManagerType::Default;
    if (!commandLineParameters.serialMode) {
        runManagerType = G4RunManagerType::MTOnly;
        cout << "Using MT run manager with " << commandLineParameters.nThreads << " threads" << endl;
    } else {
        runManagerType = G4RunManagerType::SerialOnly;
        cout << "Using serial run manager" << endl;
    }

    auto runManager = G4RunManagerFactory::CreateRunManager(runManagerType);

    if (!commandLineParameters.serialMode) {
        ROOT::EnableThreadSafety();
        runManager->SetNumberOfThreads(commandLineParameters.nThreads);
    }
#else
    cout << "Using serial run manager" << endl;
    auto runManager = new G4RunManager();
#endif

    auto detector = new DetectorConstruction(fSimulationManager);

    fSimulationManager->InitializeUserDistributions();

    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList(fSimulationManager->GetRestPhysicsLists()));

    runManager->SetUserInitialization(new ActionInitialization(fSimulationManager));

    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif

    const auto nEvents = metadata->GetNumberOfEvents();
    if (nEvents < 0) {
        cout << "Error: \"nEvents\" parameter value (" << nEvents << ") is not valid." << endl;
        exit(1);
    }

    time_t systime = time(nullptr);
    run->SetStartTimeStamp((Double_t)systime);

    cout << "Number of events: " << nEvents << endl;
    if (nEvents > 0)  // batch mode
    {
        UI->ApplyCommand("/tracking/verbose 0");
        UI->ApplyCommand("/run/initialize");
        UI->ApplyCommand("/run/beamOn " + to_string(nEvents));
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

    run->GetOutputFile()->cd();

#ifdef G4VIS_USE
    delete visManager;
#endif

    // job termination
    delete runManager;

    systime = time(nullptr);
    run->SetEndTimeStamp((Double_t)systime);
    const TString filename = TRestTools::ToAbsoluteName(run->GetOutputFileName().Data());

    // run->UpdateOutputFile(); // TODO: this line gives segfault when using MT, why?

    run->CloseFile();
    run->PrintMetadata();

    cout << "============== Generated file: " << filename << " ==============" << endl;
    auto timeEnd = chrono::steady_clock::now();
    cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(timeEnd - timeStart).count()
         << " seconds" << endl;
}
