//
// Created by lobis on 10/8/2021.
//

#include "CommandLineSetup.h"

#include <spdlog/spdlog.h>

#include <G4UImanager.hh>
#include <iostream>

using namespace std;

void CommandLineSetup::ShowUsage() {
    cout << "restG4 requires at least one parameter, the rml configuration file (-c is optional)" << endl
         << endl
         << "example: restG4 example.rml" << endl
         << endl
         << "there are other convenient optional parameters that override the ones in the rml file:" << endl
         << "\t-h or --help | show usage (this text)" << endl
         << "\t-c example.rml | specify RML file (same as calling restG4 example.rml)" << endl
         << "\t-g geometry.gdml | specify geometry file" << endl
         << "\t-i | set interactive mode (default=false)" << endl
         << "\t-s | set serial mode (no multithreading) (default=true)" << endl
         << "\t-t nThreads | set the number of threads, also enables multithreading" << endl;
}

CommandLineParameters CommandLineSetup::ProcessParameters(int argc, char** argv) {
    CommandLineParameters parameters;

    if (argc >= 2) {
        // presumably invoked as `restG4 example.rml` or `restG4 --help`
        TString argument = argv[1];
        if (argument.EqualTo("--help") | argument.EqualTo("-h")) {
            ShowUsage();
            exit(0);
        } else {
            parameters.rmlFile = argument;
        }
    }

    while (true) {
        const int option = getopt(argc, argv, "vg:c:m:o:ist:n:");
        if (option == -1) break;
        switch (option) {
            case 'c':
                if (!parameters.rmlFile.IsNull()) {
                    cout << "CommandLineParameters::ProcessParameters - Cannot specify multiple rml files "
                            " Please use at most one"
                         << endl;
                    exit(1);
                }
                parameters.rmlFile = optarg;
                break;
            case 's':
                // Serial mode
                parameters.serialMode = true;
                break;
            case 't':
                // TODO: implement this in the simulation
                // Multithreading mode
                parameters.serialMode = false;
                parameters.nThreads = std::stoi(optarg);
                if (parameters.nThreads < 1) {
                    cout << "CommandLineParameters::ProcessParameters - Number of threads must be > 0"
                         << endl;
                }
                break;
            case 'o':
                if (!parameters.outputFile.IsNull()) {
                    cout << "CommandLineParameters::ProcessParameters - Cannot specify multiple output files "
                            "via the -o flag."
                            " Please use at most one"
                         << endl;
                    exit(1);
                }
                parameters.outputFile = optarg;
                break;
            case 'i':
                parameters.interactive = true;
                break;
            case 'g':
                if (!parameters.geometryFile.IsNull()) {
                    cout << "CommandLineParameters::ProcessParameters - Cannot specify multiple geometry "
                            "files "
                            "via the -g flag. "
                            "Please use at most one"
                         << endl;
                    exit(1);
                }
                parameters.geometryFile = optarg;
                break;
            default:
                // invalid option
                cout << "Error processing command line arguments" << endl;
                exit(1);
        }
    }

    // validation
    if (parameters.rmlFile.IsNull()) {
        cout << "CommandLineParameters::ProcessParameters - Need to define an RML config file (-c)" << endl;
        ShowUsage();
        exit(1);
    }

    return parameters;
}
void CommandLineSetup::Print(const CommandLineParameters& parameters) {
    cout << "Command line parameters configuration:" << endl
         << "\t- RML file: " << parameters.rmlFile << endl
         << "\t- Execution mode: "
         << (parameters.serialMode ? "serial"
                                   : "multithreading (N = " + std::to_string(parameters.nThreads) + ")")
         << endl;
}

void CommandLineSetup::RunVisMacro() {
    spdlog::info("Setup::RunVisMacro");

    vector<string> visMacro = {
        "/control/verbose 2",
        "/run/verbose 2",
        "/run/initialize",
        "/vis/open OGL 600x600-0+0",

        "/vis/viewer/set/autoRefresh false",
        "/vis/verbose errors",
        "/vis/drawVolume",
        "/vis/viewer/set/viewpointThetaPhi 90. 0.",
        "/vis/viewer/zoom 1.5",
        //"/vis/viewer/set/style wireframe" # surface
        "/vis/scene/endOfEventAction accumulate",
        "/vis/scene/add/trajectories smooth",
        "/vis/modeling/trajectories/create/drawByCharge",
        "/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true",
        "/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2",

        "/vis/viewer/set/autoRefresh true",
        "/vis/verbose warnings",
    };

    for (const auto& command : visMacro) {
        spdlog::info("Setup::RunVisMacro - \033[1;42m{}\033[0m", command);
        G4UImanager::GetUIpointer()->ApplyCommand(command);
    }
}

void CommandLineSetup::RunGUIMacro() {
    spdlog::info("Setup::RunGUIMacro");

    vector<string> guiMacro = {
        "/gui/addMenu file File",
        "/gui/addButton file Quit exit",
        "/gui/addMenu run Run",
        R"(/gui/addButton run "beamOn 1" "/run/beamOn 1")",
        "/gui/addMenu gun Gun",
        R"(/gui/addButton gun "5 keV"   "/gun/energy 5 keV")",
        R"(/gui/addButton gun "50 keV"   "/gun/energy 50 keV")",
        R"(/gui/addButton gun "250 keV"   "/gun/energy 250 keV")",
        R"(/gui/addButton gun "1 MeV"   "/gun/energy 1 MeV")",
        // R"(/gui/addButton gun "opticalphoton"      "/gun/particle opticalphoton")",
        R"(/gui/addButton gun "e-"      "/gun/particle e-")",
        R"(/gui/addButton gun "gamma"      "/gun/particle gamma")",
        R"(/gui/addButton gun "neutron" "/gun/particle neutron")",
        R"(/gui/addButton gun "proton"  "/gun/particle proton")",
        "/gui/addMenu viewer Viewer",
        R"(/gui/addButton viewer "Set style surface" "/vis/viewer/set/style surface")",
        R"(/gui/addButton viewer "Set style wireframe" "/vis/viewer/set/style wireframe")",
        R"(/gui/addButton viewer "Refresh viewer" "/vis/viewer/refresh")",
        "/gui/addButton viewer \"Update viewer (interaction or end-of-file)\" \"/vis/viewer/update\"",
        "/gui/addButton viewer \"Flush viewer (= refresh + update)\" \"/vis/viewer/flush\"",
        R"(/gui/addButton viewer "Update scene" "/vis/scene/notifyHandlers")",
    };

    for (const auto& command : guiMacro) {
        spdlog::info("Setup::RunGUIMacro - \033[1;42m{}\033[0m", command);
        G4UImanager::GetUIpointer()->ApplyCommand(command);
    }
}
