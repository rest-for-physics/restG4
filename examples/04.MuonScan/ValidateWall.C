#include <TRestGeant4Event.h>

Int_t ValidateWall(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    if (run.GetRunTag() != "MuonsFromWall") {
        cout << "Run tag : " << run.GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'MuonsFromWall!" << endl;
        return 1;
    }

    Double_t r = 0;
    for (Int_t n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);
        Double_t x = event->GetPrimaryEventOrigin().X();
        Double_t y = event->GetPrimaryEventOrigin().Y();

        r += x * x + y * y;
    }
    r /= run.GetEntries();

    if ((Int_t)(1000. * r) / 1000 < 10000 || (Int_t)(1000. * r) / 1000 > 20000) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "R: " << (Int_t)(1000. * r) / 1000 << endl;
        return 2;
    }

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() < 350 || run.GetEntries() > 450) {
        cout << "The number of entries is not between 350 and 450!" << endl;
        cout << "Number of entries : " << run.GetEntries() << endl;
        return 3;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    // Other tests like opening other metadata classes. Detector TGeoManager, etc.

    return 0;
}
