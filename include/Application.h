
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

    int nRequestedEntries = 0;
    int timeLimitSeconds = 0;

    int runNumber = -1;

    // reference to original argc and argv necessary to pass to G4UIExecutive
    int argc;
    char** argv;
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

    void ValidateOutputFile(const std::string& outputFile) const;
};

#endif  // REST_APPLICATION_H
