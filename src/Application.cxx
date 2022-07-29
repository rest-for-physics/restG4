
#include "Application.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPRegexp.h>
#include <TROOT.h>
#include <TRestGDMLParser.h>
#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>
#include <signal.h>

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

namespace CommandLineOptions {

void ShowUsage() {
    cout << "Usage example: restG4 example.rml" << endl
         << "there are other convenient optional parameters that override the ones in the rml file:" << endl
         << "\t--help (-h) | show usage (this text)" << endl
         << "\t--config (-c) example.rml | specify RML file (same as calling restG4 example.rml)" << endl
         << "\t--output (-o) output.root | specify output file" << endl
         << "\t--events (-n) nEvents | specify number of events to be processed (overrides nEvents on rml "
            "file)"
         << endl
         << "\t--entries (-e) | specify the desired number of entries. The simulation will stop after "
            "reaching this number of saved events. Final number may be larger"
         << endl
         << "\t--time timeLimit | Sets time limit for the simulation in the format '1h20m30s', '5m20s', "
            "'30s', ... If the time limit is reached before simulation ends, it will end the simulation and "
            "save to disk all events"
         << endl
         << "\t--geometry (-g) geometry.gdml | specify geometry file" << endl
         << "\t--seed (-s) seed | specify random seed (positive integer)" << endl
         << "\t--interactive (-i) | set interactive mode (disabled by default)" << endl
         << "\t--threads (-t, -j) | set the number of threads, also enables multithreading which is disabled "
            "by default"
         << endl;
}

void PrintOptions(const Options& options) {
    cout << "Command line parameters configuration:" << endl
         << "\t- RML file: " << options.rmlFile << endl
         << (!options.outputFile.empty() ? "\t- Output file: " + options.outputFile + "\n" : "")
         << (!options.geometryFile.empty() ? "\t- Geometry file: " + options.geometryFile + "\n" : "")
         << (options.interactive ? "\t- Interactive: True\n" : "")  //
         << "\t- Execution mode: "
         << (options.nThreads == 0 ? "serial\n"
                                   : "multithreading (N = " + to_string(options.nThreads) + ")\n")
         << (options.nEvents != 0 ? "\t- Number of generated events: " + to_string(options.nEvents) + "\n"
                                  : "")
         << (options.seed != 0 ? "\t- Random seed: " + to_string(options.seed) + "\n" : "")
         << (options.nDesiredEntries != 0
                 ? "\t- Number of desired file entries: " + to_string(options.nDesiredEntries) + "\n"
                 : "")
         << (options.timeLimitSeconds != 0
                 ? "\t- Time limit: " + to_string(options.timeLimitSeconds) + " seconds\n"
                 : "")
         << endl;
}

int GetSecondsFromTimeExpression(const char* expression) {
    // expression is of the form "20h", "10m", "30s" etc.
    TPRegexp timeRegex("^(\\d+)([hms])$$");
    TObjArray* subStrL = timeRegex.MatchS(expression);
    const Int_t nrSubStr = subStrL->GetLast() + 1;
    if (nrSubStr > 2) {
        const int time = stoi(((TObjString*)subStrL->At(1))->GetString().Data());
        const TString modifier = ((TObjString*)subStrL->At(2))->GetString();

        if (modifier == "h") {
            return time * 60 * 60;
        } else if (modifier == "m") {
            return time * 60;
        } else if (modifier == "s") {
            return time;
        }
    }

    return 0;
}

int GetSecondsFromFullTimeExpression(const char* expression) {
    // expression is of the form "1h20m30s", "1h", "20m30s", "10s" etc.
    int seconds = 0;
    TPRegexp fullTimeRegex("^(\\d+h)?(\\d+m)?(\\d+s)?$");
    TObjArray* subStrL = fullTimeRegex.MatchS(expression);

    for (int i = 0; i < 3; i++) {
        auto obj = (TObjString*)subStrL->At(i + 1);
        if (obj != nullptr) {
            seconds += GetSecondsFromTimeExpression(obj->GetString().Data());
        }
    }

    return seconds;
}

Options ProcessCommandLineOptions(int argc, char* const argv[]) {
    Options options;
    options.argc = argc;
    options.argv = const_cast<char**>(argv);

    if (argc < 2) {
        // Invoked without parameter
        ShowUsage();
        exit(0);
    }

    // See https://cplusplus.com/articles/DEN36Up4/
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            ShowUsage();
            exit(0);
        } else if ((arg == "-c") || (arg == "--config")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.rmlFile =
                    argv[++i];  // Increment 'i' so we don't get the argument as the next argv[i].
            } else {
                cerr << "--config option requires one argument" << endl;
                exit(1);
            }
        } else if ((arg == "-o") || (arg == "--output")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.outputFile =
                    argv[++i];  // Increment 'i' so we don't get the argument as the next argv[i].
            } else {
                cerr << "--output option requires one argument" << endl;
                exit(1);
            }
        } else if ((arg == "-g") || (arg == "--geometry")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.geometryFile =
                    argv[++i];  // Increment 'i' so we don't get the argument as the next argv[i].
            } else {
                cerr << "--geometry option requires one argument" << endl;
                exit(1);
            }
        } else if ((arg == "-i") || (arg == "--interactive")) {
            options.interactive = true;
            // TODO: not yet implemented
            cout << "--interactive option not yet implemented" << endl;
            exit(1);
        } else if ((arg == "-j") || (arg == "--threads") || (arg == "-t")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.nThreads =
                    stoi(argv[++i]);  // Increment 'i' so we don't get the argument as the next argv[i].
                if (options.nThreads < 1) {
                    cout << "--threads option error: number of threads must be > 0" << endl;
                    exit(1);
                }
            } else {
                cerr << "--threads option requires one argument." << endl;
                exit(1);
            }
        } else if ((arg == "-n") || (arg == "--events")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.nEvents =
                    stoi(argv[++i]);  // Increment 'i' so we don't get the argument as the next argv[i].
                if (options.nEvents <= 0) {
                    cout << "--events option error: number of events must be > 0" << endl;
                    exit(1);
                }
            } else {
                cerr << "--events option requires one argument." << endl;
                exit(1);
            }
        } else if ((arg == "-e") || (arg == "--entries")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.nDesiredEntries =
                    stoi(argv[++i]);  // Increment 'i' so we don't get the argument as the next argv[i].
                if (options.nDesiredEntries <= 0) {
                    cout << "--entries option error: number of entries must be > 0" << endl;
                    exit(1);
                }
            } else {
                cerr << "--entries option requires one argument." << endl;
                exit(1);
            }
        } else if ((arg == "-s") || (arg == "--seed")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.seed =
                    stoi(argv[++i]);  // Increment 'i' so we don't get the argument as the next argv[i].
                if (options.seed <= 0) {
                    cout << "--seed option error: seed must be positive number" << endl;
                    exit(1);
                }
            } else {
                cerr << "--seed option requires one argument." << endl;
                exit(1);
            }
        } else if (arg == "--time") {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                options.timeLimitSeconds = GetSecondsFromFullTimeExpression(
                    argv[++i]);  // Increment 'i' so we don't get the argument as the next argv[i].
                if (options.timeLimitSeconds <= 0) {
                    cout << "--time option error: time limit must be of the format 1h20m30s, 10m20s, 1h, etc."
                         << endl;
                    exit(1);
                }
            } else {
                cerr << "--time option requires one argument." << endl;
                exit(1);
            }
        } else {
            options.rmlFile = argv[i];  // invoked as restG4 <rmlFile>, restG4 -j 4 <rmlFile>, etc.
        }
    }

    if (options.rmlFile.empty()) {
        cerr << "Input RML file not specified" << endl;
        exit(1);
    }

    return options;
}

}  // namespace CommandLineOptions

int interruptSignalHandler(const int, void* ptr) {
    // See https://stackoverflow.com/a/43400143/11776908
    cout << "Stopping Run! Program was manually stopped by user (CTRL+C)!" << endl;
    const auto manager = (SimulationManager*)(ptr);
    manager->StopSimulation();
    return 0;
}

void Application::Run(const CommandLineOptions::Options& options) {
    signal(SIGINT, (void (*)(int))interruptSignalHandler);

    const auto originalDirectory = filesystem::current_path();

    cout << "Current working directory: " << originalDirectory << endl;

    CommandLineOptions::PrintOptions(options);

    // Separating relative path and pure RML filename
    const char* inputConfigFile = options.rmlFile.c_str();

    if (!TRestTools::CheckFileIsAccessible(inputConfigFile)) {
        cerr << "Input RML file " << filesystem::weakly_canonical(inputConfigFile)
             << " not found, please check file name!" << endl;
        exit(1);
    }

    const auto [inputRmlPath, inputRmlClean] = TRestTools::SeparatePathAndName(inputConfigFile);

    if (!filesystem::path(inputRmlPath).empty()) {
        filesystem::current_path(inputRmlPath);
    }

    auto metadata = new TRestGeant4Metadata(inputRmlClean.c_str());
    fSimulationManager.SetRestMetadata(metadata);

    metadata->SetGeant4Version(TRestTools::Execute("geant4-config --version"));

    if (options.seed != 0) {
        metadata->SetSeed(options.seed);
    }
    constexpr auto maxPrimariesAllowed = 2147483647;
    if (options.nEvents != 0) {
        metadata->SetNumberOfEvents(options.nEvents);
    }
    if (metadata->GetNumberOfEvents() == 0) {
        metadata->SetNumberOfEvents(maxPrimariesAllowed);
    }
    if (options.nDesiredEntries != 0) {
        metadata->SetNumberOfDesiredEntries(options.nDesiredEntries);
    }
    if (options.timeLimitSeconds != 0) {
        metadata->SetSimulationMaxTimeSeconds(options.timeLimitSeconds);
    }
    if (!options.geometryFile.empty()) {
        metadata->SetGdmlFilename(options.geometryFile);
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

    metadata->PrintMetadata();

    auto physicsLists = new TRestGeant4PhysicsLists(inputRmlClean.c_str());
    fSimulationManager.SetRestPhysicsLists(physicsLists);

    auto run = new TRestRun();
    fSimulationManager.SetRestRun(run);

    run->LoadConfigFromFile(inputRmlClean);

    if (!options.outputFile.empty()) {
        run->SetOutputFileName(options.outputFile);
    }

    filesystem::current_path(originalDirectory);

    TString runTag = run->GetRunTag();
    if (runTag == "Null" || runTag == "") {
        run->SetRunTag(metadata->GetTitle());
    }

    run->SetRunType("restG4");

    run->AddMetadata(fSimulationManager.GetRestMetadata());
    run->AddMetadata(fSimulationManager.GetRestPhysicsLists());

    run->PrintMetadata();

    run->FormOutputFile();
    run->GetOutputFile()->cd();

    run->AddEventBranch(&fSimulationManager.fEvent);

    // choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    long seed = metadata->GetSeed();
    CLHEP::HepRandom::setTheSeed(seed);

    G4VSteppingVerbose::SetInstance(new SteppingVerbose(&fSimulationManager));

#ifndef GEANT4_WITHOUT_G4RunManagerFactory
    auto runManagerType = G4RunManagerType::Default;
    const bool serialMode = options.nThreads == 0;
    if (serialMode) {
        runManagerType = G4RunManagerType::SerialOnly;
        cout << "Using serial run manager" << endl;
    } else {
        runManagerType = G4RunManagerType::MTOnly;
        cout << "Using MT run manager with " << options.nThreads << " threads" << endl;
    }

    auto runManager = G4RunManagerFactory::CreateRunManager(runManagerType);

    if (!serialMode) {
        ROOT::EnableThreadSafety();
        runManager->SetNumberOfThreads(options.nThreads);
    }
#else
    cout << "Using serial run manager" << endl;
    auto runManager = new G4RunManager();
#endif

    auto detector = new DetectorConstruction(&fSimulationManager);

    fSimulationManager.InitializeUserDistributions();

    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList(fSimulationManager.GetRestPhysicsLists()));
    runManager->SetUserInitialization(new ActionInitialization(&fSimulationManager));

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
        auto ui = new G4UIExecutive(options.argc, options.argv);
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

    const auto nEntries = run->GetEntries();

    run->UpdateOutputFile();
    run->CloseFile();

    auto geometry = gdml->CreateGeoManager();
    WriteGeometry(geometry, run->GetOutputFileName());
    delete geometry;

    metadata->PrintMetadata();
    run->PrintMetadata();

    cout << "\t- Total simulation time is " << fSimulationManager.GetElapsedTime() << " seconds, " << nEvents
         << " processed events (" << nEvents / fSimulationManager.GetElapsedTime() << " per second) and "
         << nEntries << " events saved to output file (" << nEntries / fSimulationManager.GetElapsedTime()
         << " per second)" << endl;
    cout << "\t- Output file: " << filename << endl << endl;
}

void Application::WriteGeometry(TGeoManager* geometry, const char* filename, const char* option) {
    auto file = TFile::Open(filename, option);
    file->cd();
    cout << "Application::WriteGeometry - Writing geometry into '" << filename << "'" << endl;
    if (!geometry) {
        cout << "Application::WriteGeometry - Error - Unable to write geometry into file" << endl;
        exit(1);
    }
    geometry->Write("Geometry");

    file->Close();
}
