
#ifndef REST_APPLICATION_H
#define REST_APPLICATION_H

#include "SimulationManager.h"

class G4VisManager;
class G4UIExecutive;
class G4RunManager;

class TGeoManager;

namespace CommandLineOptions {
struct Options {
    std::string rmlFile{};
    std::string outputFile{};
    std::string geometryFile{};

    bool interactive = false;

    int nThreads = 0;

    int nEvents = 0;
    Long_t seed = 0;

    int nDesiredEntries = 0;
    int timeLimitSeconds = 0;
};

Options ProcessCommandLineOptions(int argc, char* const argv[]);
void PrintOptions(const Options& options);
void ShowUsage();

}  // namespace CommandLineOptions

class Application {
   public:
    void Run(const CommandLineOptions::Options& options);

    ~Application() = default;

   private:
    SimulationManager fSimulationManager;

    void WriteGeometry(TGeoManager* geometry, const char* filename, const char* option = "UPDATE");
};

#endif  // REST_APPLICATION_H
