
#include "Application.h"

using namespace std;

constexpr auto banner =
    "\n"
    "                         888     .d8888b.      d8888 \n"
    "                         888    d88P  Y88b    d8P888 \n"
    "                         888    888    888   d8P 888 \n"
    "888d888 .d88b.  .d8888b  888888 888         d8P  888 \n"
    "888P\"  d8P  Y8b 88K      888    888  88888 d88   888 \n"
    "888    88888888 \"Y8888b. 888    888    888 8888888888\n"
    "888    Y8b.          X88 Y88b.  Y88b  d88P       888 \n"
    "888     \"Y8888   88888P'  \"Y888  \"Y8888P88       888 \n"
    "\n"
    "🔍 Source code: https://github.com/rest-for-physics/restG4\n"
    "🐞 Bug reports: https://github.com/rest-for-physics/restG4/issues\n";

int main(int argc, char* argv[]) {
    cout << banner << endl;

    Application app;

    CommandLineOptions::Options options = CommandLineOptions::ProcessCommandLineOptions(argc, argv);

    app.Run(options);
}
