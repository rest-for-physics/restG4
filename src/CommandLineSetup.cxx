//
// Created by lobis on 10/8/2021.
//

#include "CommandLineSetup.h"

#include <getopt.h>
#include <unistd.h>

#include <iostream>

using namespace std;

void CommandLineSetup::ShowUsage() {
    cout << "restG4 requires at least one parameter, the rml configuration file" << endl
         << endl
         << "example: restG4 example.rml" << endl
         << endl
         << "there are other convenient optional parameters that override the ones in the rml file:" << endl
         << "\t-h or --help | show usage (this text)" << endl
         << "\t-c example.rml | specify RML file (same as calling restG4 example.rml)" << endl
         << "\t-o output.root | specify output file" << endl
         << "\t-g geometry.gdml | specify geometry file" << endl
         << "\t-i | set interactive mode (default=false)" << endl
         << "\t-s | set serial mode (no multithreading) (default=true)" << endl
         << "\t-t nThreads | set the number of threads, also enables multithreading" << endl;
}

CommandLineParameters CommandLineSetup::ProcessParameters(int argc, char** argv) {
    CommandLineParameters parameters;

    parameters.argc = argc;
    parameters.argv = argv;

    if (argc >= 2) {
        // presumably invoked as `restG4 example.rml` or `restG4 --help`
        TString argument = argv[1];
        if (argument.EqualTo("--help") || argument.EqualTo("-h")) {
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
                parameters.nThreads = stoi(optarg);
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
                            "files via the -g flag. Please use at most one"
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
         << (!parameters.outputFile.IsNull() ? "\t- Output file: " + parameters.outputFile + "\n" : "")
         << (!parameters.geometryFile.IsNull() ? "\t- Geometry file: " + parameters.geometryFile + "\n" : "")
         << (parameters.interactive ? "\t- Interactive: True\n" : "")  //
         << "\t- Execution mode: "
         << (parameters.serialMode ? "serial" : "multithreading (N = " + to_string(parameters.nThreads) + ")")
         << endl;
}
