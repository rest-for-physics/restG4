//
// Created by lobis on 10/8/2021.
//

#include "CommandLineSetup.h"

#include <TObjArray.h>
#include <TObjString.h>
#include <TPRegexp.h>
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
         << "\t-n nEvents | specify number of events to be generated (overrides nEvents on rml file). "
            "Incompatible with '-N' option"
         << endl
         << "\t-N nEventsOnFile | specify desired number of entries after simulation is completed. It will "
            "adjust 'nEvents' accordingly to get this number. Final number of entries may be larger. "
            "Incompatible with '-n' option"
         << endl
         << "\t-S timeLimit | Sets time limit for the simulation. If the time limit is reached before "
            "simulation ends, it will end the simulation and save to disk all events. Follows expression "
            "like '1h20m30s', '5m20s', '30s', ..."
         << endl
         << "\t-g geometry.gdml | specify geometry file" << endl
         << "\t-i | set interactive mode (default=false)" << endl
         << "\t-s | set serial mode (no multithreading) (default=true)" << endl
         << "\t-t nThreads | set the number of threads, also enables multithreading" << endl;
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

CommandLineParameters CommandLineSetup::ProcessParameters(int argc, char** argv) {
    CommandLineParameters parameters;

    parameters.cmdArgc = argc;
    parameters.cmdArgv = argv;

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
        const int option = getopt(argc, argv, "vg:c:m:o:ist:n:N:S:");
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
                parameters.nThreads = 0;
                break;
            case 't':
                // TODO: implement this in the simulation
                // Multithreading mode
                parameters.nThreads = stoi(optarg);
                if (parameters.nThreads < 1) {
                    cout << "CommandLineParameters::ProcessParameters - Number of threads must be > 0"
                         << endl;
                    exit(1);
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
            case 'n':
                if (parameters.nDesiredEntries != 0) {
                    cout << "CommandLineParameters::ProcessParameters - '-n' and '-N' options are mutually "
                            "incompatible"
                         << endl;
                    exit(1);
                }
                parameters.nEvents = stoi(optarg);
                if (parameters.nEvents < 1) {
                    cout
                        << "CommandLineParameters::ProcessParameters - number of generated events must be > 0"
                        << endl;
                    exit(1);
                }
                break;
            case 'N':
                if (parameters.nEvents != 0) {
                    cout << "CommandLineParameters::ProcessParameters - '-n' and '-N' options are mutually "
                            "incompatible"
                         << endl;
                    exit(1);
                }
                parameters.nDesiredEntries = stoi(optarg);
                if (parameters.nDesiredEntries < 1) {
                    cout << "CommandLineParameters::ProcessParameters - number of desired entries must be > 0"
                         << endl;
                    exit(1);
                }
                break;
            case 'S':
                parameters.timeLimitSeconds = GetSecondsFromFullTimeExpression(optarg);
                if (parameters.timeLimitSeconds == 0) {
                    cout << "CommandLineParameters::ProcessParameters - time regex failed to match for "
                            "argument '"
                         << optarg << "'" << endl;
                    exit(1);
                }
                break;
            case 'i':
                parameters.interactive = true;
                break;
            case 'g':
                // TODO: implement
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
         << (parameters.nThreads == 0 ? "serial\n"
                                      : "multithreading (N = " + to_string(parameters.nThreads) + ")\n")
         << (parameters.nEvents != 0
                 ? "\t- Number of generated events: " + to_string(parameters.nEvents) + "\n"
                 : "")
         << (parameters.nDesiredEntries != 0
                 ? "\t- Number of desired file entries: " + to_string(parameters.nDesiredEntries) + "\n"
                 : "")
         << (parameters.timeLimitSeconds != 0
                 ? "\t- Time limit: " + to_string(parameters.timeLimitSeconds) + " seconds\n"
                 : "")
         << endl;
}
