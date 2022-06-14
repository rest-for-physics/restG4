
#include <Application.h>
#include <TRestGeant4Metadata.h>
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