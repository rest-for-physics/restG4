
#include <Application.h>
#include <CommandLineSetup.h>
#include <TROOT.h>
#include <gtest/gtest.h>

#include <filesystem>

namespace fs = std::filesystem;

using namespace std;

const auto examplesPath = fs::path(__FILE__).parent_path().parent_path().parent_path() / "examples";

TEST(restG4, CheckExampleFiles) {
    cout << "Examples files path: " << examplesPath << endl;

    // Check dir exists and is a non-empty directory
    EXPECT_TRUE(fs::is_directory(examplesPath));
    EXPECT_TRUE(!fs::is_empty(examplesPath));
}

TEST(restG4, Example_01_NLDBD) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "01.NLDBD";
    fs::current_path(thisExamplePath);

    CommandLineParameters parameters;
    parameters.rmlFile = "NLDBD.rml";
    parameters.outputFile =
        thisExamplePath / "NLDBD_simulation.root";  // TODO: fix not working with local path

    Application app;
    app.Run(parameters);

    // Run validation macro
    const TString macro(thisExamplePath / "Validate.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("Validate(\"%s\")", parameters.outputFile.Data()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}
