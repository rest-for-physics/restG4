
#include "Application.h"

using namespace std;

int main(int argc, char* argv[]) {
    Application app;

    CommandLineOptions::Options options = CommandLineOptions::ProcessCommandLineOptions(argc, argv);

    app.Run(options);
}
