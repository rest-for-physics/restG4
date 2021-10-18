//
// Created by lobis on 10/8/2021.
//

#ifndef REST_COMMANDLINESETUP_H
#define REST_COMMANDLINESETUP_H

#include <TString.h>

class CommandLineParameters {
   public:
    TString rmlFile;           // -c (also default argument, does not need '-c') REQUIRED
    TString outputFile;        // -o output file OPTIONAL
    TString geometryFile;      // -g UNIQUE OPTIONAL
    bool interactive = false;  // -i FLAG (NO ARGUMENT) DEFAULT=FALSE OPTIONAL
    bool serialMode = true;    // -s FLAG (NO ARGUMENT) DEFAULT=TRUE OPTIONAL
    int nThreads = 0;          // -t Number of threads, only active if serialMode is set to false OPTIONAL

    CommandLineParameters() = default;
};

class CommandLineSetup {
   public:
    CommandLineSetup() = delete;

    static CommandLineParameters ProcessParameters(int argc, char* argv[]);

    static void ShowUsage();
    static void Print(const CommandLineParameters&);
};

#endif  // REST_COMMANDLINESETUP_H