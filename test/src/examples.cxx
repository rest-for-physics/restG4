
#include <Application.h>
#include <CommandLineSetup.h>
#include <TROOT.h>
#include <TRestRun.h>
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

    {  // Run simulation
        Application app;
        app.Run(parameters);
    }

    {  // Run validation macro
        const TString macro(thisExamplePath / "Validate.C");
        gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
        int error = 0;
        const int result =
            gROOT->ProcessLine(TString::Format("Validate(\"%s\")", parameters.outputFile.Data()), &error);
        EXPECT_EQ(error, 0);
        EXPECT_EQ(result, 0);
    }

    fs::current_path(originalPath);
}

TEST(restG4, TRestGeant4GeometryInfo_TRestGeant4PhysicsInfo) {
    // Test "TRestGeant4GeometryInfo" and "TRestGeant4PhysicsInfo" even though its from Geant4Lib, we need a
    // simulation file, so we placed the test here

    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "01.NLDBD";
    const auto resultsFile = thisExamplePath / "NLDBD_simulation.root";

    cout << "results file: " << resultsFile << endl;

    // If previous example was ran we can use generated file, otherwise we generate it again
    if (!TRestTools::CheckFileIsAccessible(resultsFile)) {
        cout << "Results file not found, generating file..." << endl;
        fs::current_path(thisExamplePath);

        CommandLineParameters parameters;
        parameters.rmlFile = "NLDBD.rml";
        parameters.outputFile = resultsFile;

        Application app;
        app.Run(parameters);

        fs::current_path(originalPath);
    } else {
        cout << "Results file found, skipping generation" << endl;
    }

    TRestRun run(resultsFile);
    // Test `TRestGeant4Metadata::GetUnambiguousGlobalInstance`
    auto geant4Metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    EXPECT_EQ(geant4Metadata != nullptr, true);

    const auto& geometryInfo = geant4Metadata->GetGeant4GeometryInfo();

    cout << "Printing geometry info" << endl;
    geometryInfo.Print();

    EXPECT_EQ(geometryInfo.GetAllLogicalVolumes().size() == 4, true);
    EXPECT_EQ(geometryInfo.GetAllPhysicalVolumes().size() == 4, true);

    for (const auto& logicalVolume : geometryInfo.GetAllLogicalVolumes()) {
        if (logicalVolume == "World") {
            const auto& physicals = geometryInfo.GetAllPhysicalVolumesFromLogical(logicalVolume);
            EXPECT_EQ(physicals.size() == 1, true);
            const auto& physical = physicals[0];
            EXPECT_EQ(physical == "World_PV", true);
        } else if (logicalVolume == "gasVolume") {
            const auto& physicals = geometryInfo.GetAllPhysicalVolumesFromLogical(logicalVolume);
            EXPECT_EQ(physicals.size() == 1, true);
            const auto& physical = physicals[0];
            EXPECT_EQ(physical == "gas", true);
        } else if (logicalVolume == "vesselVolume") {
            const auto& physicals = geometryInfo.GetAllPhysicalVolumesFromLogical(logicalVolume);
            EXPECT_EQ(physicals.size() == 1, true);
            const auto& physical = physicals[0];
            EXPECT_EQ(physical == "vessel", true);
        } else if (logicalVolume == "waterTankVolume") {
            const auto& physicals = geometryInfo.GetAllPhysicalVolumesFromLogical(logicalVolume);
            EXPECT_EQ(physicals.size() == 1, true);
            const auto& physical = physicals[0];
            EXPECT_EQ(physical == "waterTank", true);
        } else {
            cout << "unexpected logical volume found!" << endl;
            exit(1);
        }
    }

    const auto& physicsInfo = geant4Metadata->GetGeant4PhysicsInfo();

    cout << "Printing physics info" << endl;
    physicsInfo.Print();

    const auto particles = physicsInfo.GetAllParticles();
    EXPECT_EQ(particles.size() == 2, true);
    EXPECT_EQ(particles.count("e-") > 0, true);
    EXPECT_EQ(particles.count("gamma") > 0, true);

    const auto processes = physicsInfo.GetAllProcesses();
    EXPECT_EQ(processes.count("Init") > 0, true);
    EXPECT_EQ(processes.count("Transportation") > 0, true);
    EXPECT_EQ(processes.count("compt") > 0, true);
}

TEST(restG4, Example_04_Muons) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "04.MuonScan";
    fs::current_path(thisExamplePath);

    CommandLineParameters parameters;
    parameters.rmlFile = "CosmicMuonsFromWall.rml";
    parameters.outputFile = thisExamplePath / "muons.root";  // TODO: fix not working with local path

    Application app;
    app.Run(parameters);

    // Run validation macro
    const TString macro(thisExamplePath / "ValidateWall.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("ValidateWall(\"%s\")", parameters.outputFile.Data()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);

    // use output file to check additional things
    TRestRun run(parameters.outputFile.Data());
    cout << "Number of entries: " << run.GetEntries() << endl;
}

TEST(restG4, Example_05_PandaX) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "05.PandaXIII";
    fs::current_path(thisExamplePath);

    CommandLineParameters parameters;
    parameters.rmlFile = "Xe136bb0n.rml";
    parameters.outputFile = thisExamplePath / "Xe136bb0n.root";  // TODO: fix not working with local path

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

TEST(restG4, Example_06_IonRecoils) {
    GTEST_SKIP_("Not working yet, needs fixing");
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "06.IonRecoils";
    fs::current_path(thisExamplePath);

    CommandLineParameters parameters;
    parameters.rmlFile = "recoils.rml";
    parameters.outputFile = thisExamplePath / "recoils.root";

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

TEST(restG4, Example_07_Decay_FullChain) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "07.FullChainDecay";
    fs::current_path(thisExamplePath);

    CommandLineParameters parameters;
    parameters.rmlFile = "fullChain.rml";
    parameters.outputFile = thisExamplePath / "fullChain.root";  // TODO: fix not working with local path

    Application app;
    app.Run(parameters);

    // print processes
    TRestRun run(parameters.outputFile.Data());
    auto geant4Metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    EXPECT_EQ(geant4Metadata != nullptr, true);

    const auto& physicsInfo = geant4Metadata->GetGeant4PhysicsInfo();

    cout << "Printing physics info" << endl;
    physicsInfo.Print();

    // Run validation macro
    const TString macro(thisExamplePath / "Validate.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("Validate(\"%s\", %d)", parameters.outputFile.Data(), 16), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}
