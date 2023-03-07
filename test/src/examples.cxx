
#include <Application.h>
#include <TGeoManager.h>
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

    CommandLineOptions::Options options;
    options.rmlFile = "NLDBD.rml";
    options.outputFile = thisExamplePath / "NLDBD_simulation.root";  // TODO: fix not working with local path

    {  // Run simulation
        Application app;
        app.Run(options);
    }

    {  // Run validation macro
        const TString macro(thisExamplePath / "Validate.C");
        gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
        int error = 0;
        const int result =
            gROOT->ProcessLine(TString::Format("Validate(\"%s\")", options.outputFile.c_str()), &error);
        EXPECT_EQ(error, 0);
        EXPECT_EQ(result, 0);
    }

    fs::current_path(originalPath);
}

TEST(restG4, Metadata) {
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

        CommandLineOptions::Options options;
        options.rmlFile = "NLDBD.rml";
        options.outputFile = resultsFile;

        Application app;
        app.Run(options);

        fs::current_path(originalPath);
    } else {
        cout << "Results file found, skipping generation" << endl;
    }

    TRestRun run(resultsFile);

    /* Check TGeoManager is present on file */
    const TGeoManager* geometry = run.GetInputFile()->Get<TGeoManager>("Geometry");
    EXPECT_EQ(geometry != nullptr, true);
    delete geometry;

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

    CommandLineOptions::Options options;
    options.rmlFile = "CosmicMuonsFromWall.rml";
    options.outputFile = thisExamplePath / "CosmicMuonsFromWall.root";

    Application app;
    app.Run(options);

    // Run validation macro
    const TString macro(thisExamplePath / "ValidateCosmicMuonsFromWall.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result = gROOT->ProcessLine(
        TString::Format("ValidateCosmicMuonsFromWall(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);

    // use output file to check additional things
    TRestRun run(options.outputFile.c_str());
    cout << "Number of entries: " << run.GetEntries() << endl;
}

TEST(restG4, MergeFiles) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "04.MuonScan";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "CosmicMuonsFromWall.rml";

    // create "merge" directory
    auto mergeDirectory = fs::create_directory(thisExamplePath / "merge");
    options.nEvents = 1000;

    options.outputFile = mergeDirectory / "muons1.root";
    {
        Application app;
        app.Run(options);
    }

    options.nEvents = 500;
    options.outputFile = mergeDirectory / "muons2.root";
    {
        Application app;
        app.Run(options);
    }

    options.nEvents = 200;
    options.outputFile = mergeDirectory / "muons3.root";
    {
        Application app;
        app.Run(options);
    }

    // run system command to merge files
    const auto command = "restMergeFiles merge.root " + mergeDirectory.string();
    cout << "Running command: " << command << endl;
    const auto result = system(command.c_str());
    EXPECT_EQ(result, 0);

    int nEntries = 0;
    {
        TRestRun run("muons1.root");
        nEntries += run.GetEntries();
    }
    {
        TRestRun run("muons2.root");
        nEntries += run.GetEntries();
    }
    {
        TRestRun run("muons3.root");
        nEntries += run.GetEntries();
    }

    cout << "Computed number of entries: " << nEntries << endl;

    TRestRun run("merge.root");
    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    cout << "Number of entries: " << run.GetEntries() << endl;
    cout << "Number of primaries: " << metadata->GetNumberOfEvents() << endl;

    EXPECT_EQ(metadata->GetNumberOfEvents(), 1700);
    EXPECT_EQ(run.GetEntries(), nEntries);
}

TEST(restG4, Example_04_Muons_MT) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "04.MuonScan";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "CosmicMuonsFromWall.rml";
    options.outputFile = thisExamplePath / "CosmicMuonsFromWall.root";

    options.nThreads = 4;

    Application app;
    app.Run(options);

    // Run validation macro
    const TString macro(thisExamplePath / "ValidateCosmicMuonsFromWall.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result = gROOT->ProcessLine(
        TString::Format("ValidateCosmicMuonsFromWall(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);

    // use output file to check additional things
    TRestRun run(options.outputFile.c_str());
    cout << "Number of entries: " << run.GetEntries() << endl;
}

TEST(restG4, Example_05_PandaX) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "05.PandaXIII";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "Xe136bb0n.rml";
    options.outputFile = thisExamplePath / "Xe136bb0n.root";  // TODO: fix not working with local path

    Application app;
    app.Run(options);

    // Run validation macro
    const TString macro(thisExamplePath / "Validate.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("Validate(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_06_IonRecoils) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "06.IonRecoils";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "recoils.rml";
    options.outputFile = thisExamplePath / "recoils.root";

    Application app;
    app.Run(options);

    // Run validation macro
    const TString macro(thisExamplePath / "Validate.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("Validate(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_07_Decay_FullChain) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "07.FullChainDecay";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "fullChain.rml";
    options.outputFile = thisExamplePath / "fullChain.root";  // TODO: fix not working with local path

    Application app;
    app.Run(options);

    // print processes
    TRestRun run(options.outputFile.c_str());
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
        gROOT->ProcessLine(TString::Format("Validate(\"%s\", %d)", options.outputFile.c_str(), 15), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_09_Pb210_Shield) {
    GTEST_SKIP_("This test should work, but we skip it because it takes too long");
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "09.Pb210_Shield";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "Pb210.rml";
    options.outputFile = thisExamplePath / "shielding.root";  // TODO: fix not working with local path

    Application app;
    app.Run(options);

    TRestRun run(options.outputFile.c_str());

    // Run validation macro
    const TString macro(thisExamplePath / "Validate.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("Validate(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_10_Geometry) {
    // cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "10.Geometries";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "Assembly.rml";
    options.outputFile = thisExamplePath / "geometries.root";

    Application app;
    app.Run(options);

    TRestRun run(options.outputFile);

    // Run validation macro
    const TString macro(thisExamplePath / "Validate.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("Validate(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_12_Generators_CosineSquared) {
    //  cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "12.Generators";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "CosineSquaredCircle.rml";
    options.outputFile = thisExamplePath / "CosineSquaredCircle.root";

    Application app;
    app.Run(options);

    TRestRun run(options.outputFile);

    // Run validation macro
    const TString macro(thisExamplePath / "ValidateCosineSquared.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result = gROOT->ProcessLine(
        TString::Format("ValidateCosineSquared(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_12_Generators_EnergyAndAngularRange) {
    //  cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "12.Generators";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "EnergyAndAngularRange.rml";
    options.outputFile = thisExamplePath / "EnergyAndAngularRange.root";

    Application app;
    app.Run(options);

    TRestRun run(options.outputFile);

    // Run validation macro
    const TString macro(thisExamplePath / "ValidateEnergyAndAngularRange.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result = gROOT->ProcessLine(
        TString::Format("ValidateEnergyAndAngularRange(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_12_Generators_EnergyAndAngularCorrelated) {
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "12.Generators";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "CosmicMuonsEnergyAngularCorrelated.rml";
    options.outputFile = thisExamplePath / "CosmicMuonsEnergyAngularCorrelated.root";

    Application app;
    app.Run(options);

    // Run validation macro
    const TString macro(thisExamplePath / "ValidateCosmicMuonsEnergyAngularCorrelated.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result = gROOT->ProcessLine(
        TString::Format("ValidateCosmicMuonsEnergyAngularCorrelated(\"%s\")", options.outputFile.c_str()),
        &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}

TEST(restG4, Example_13_IAXO_Neutrons) {
    //  cd into example
    const auto originalPath = fs::current_path();
    const auto thisExamplePath = examplesPath / "13.IAXO";
    fs::current_path(thisExamplePath);

    CommandLineOptions::Options options;
    options.rmlFile = "Neutrons.rml";
    options.outputFile = thisExamplePath / "Neutrons.root";
    options.nRequestedEntries = 1;

    Application app;
    app.Run(options);

    TRestRun run(options.outputFile);

    // Run validation macro
    const TString macro(thisExamplePath / "Validate.C");
    gROOT->ProcessLine(TString::Format(".L %s", macro.Data()));  // Load macro
    int error = 0;
    const int result =
        gROOT->ProcessLine(TString::Format("Validate(\"%s\")", options.outputFile.c_str()), &error);
    EXPECT_EQ(error, 0);
    EXPECT_EQ(result, 0);

    fs::current_path(originalPath);
}
