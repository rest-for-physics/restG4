
#include "Application.h"

using namespace std;

constexpr auto banner =
    "                             #####  #       \n"
    " #####  ######  ####  ##### #     # #    #  \n"
    " #    # #      #        #   #       #    #  \n"
    " #    # #####   ####    #   #  #### #    #  \n"
    " #####  #           #   #   #     # ####### \n"
    " #   #  #      #    #   #   #     #      #  \n"
    " #    # ######  ####    #    #####       #  \n"
    "\n"
    "🔍 Source code: https://github.com/rest-for-physics/restG4\n"
    "🐞 Bug reports: https://github.com/rest-for-physics/restG4/issues\n";

int main(int argc, char* argv[]) {
    cout << banner << endl;

    Application app;

    CommandLineOptions::Options options = CommandLineOptions::ProcessCommandLineOptions(argc, argv);

    app.Run(options);
}
