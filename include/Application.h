
#ifndef REST_APPLICATION_H
#define REST_APPLICATION_H

#include "CommandLineSetup.h"
#include "SimulationManager.h"

class G4VisManager;
class G4UIExecutive;
class G4RunManager;

class TGeoManager;

class Application {
   public:
    inline static CommandLineParameters ProcessCommandLineParameters(int argc, char** argv) {
        return CommandLineSetup::ProcessParameters(argc, argv);
    }

    void Run(const CommandLineParameters& parameters);
    inline void Run(int argc, char** argv) { Run(ProcessCommandLineParameters(argc, argv)); }

    ~Application() = default;

   private:
    SimulationManager fSimulationManager;

    void WriteGeometry(TGeoManager* geometry, const char* filename, const char* option = "UPDATE");
};

#endif  // REST_APPLICATION_H
