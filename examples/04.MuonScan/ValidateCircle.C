#include <TRestGeant4Event.h>

Int_t ValidateCircle(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    if (run.GetRunTag() != "MuonsFromCircle") {
        cout << "Run tag: " << run.GetRunTag() << endl;
        cout << "The run tag of the basic validation test should be 'MuonsFromCircle'" << endl;
        return 1;
    }

    double radiusAverage = 0, radiusMin = TMath::Infinity(), radiusMax = 0;
    constexpr double radiusAverageRef = 2.92, radiusMinRef = 0.0, radiusMaxRef = 4.0;
    constexpr double tolerance = 0.1;

    for (Int_t n = 0; n < run.GetEntries(); n++) {
        run.GetEntry(n);
        Double_t x = event->GetPrimaryEventOrigin().X();
        Double_t z = event->GetPrimaryEventOrigin().Z();
        const auto r = TMath::Sqrt(x * x + z * z) / 100;  // cm
        radiusAverage += r / run.GetEntries();
        if (r < radiusMin) {
            radiusMin = r;
        }
        if (r > radiusMax) {
            radiusMax = r;
        }
    }

    cout << "Average radius (cm): " << radiusAverage << endl;
    cout << "Minimum radius (cm): " << radiusMin << endl;
    cout << "Maximum radius (cm): " << radiusMax << endl;

    if (TMath::Abs(radiusAverage - radiusAverageRef) > tolerance) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "radiusAverage (cm): " << radiusAverage << endl;
        return 3;
    }
    if (TMath::Abs(radiusMin - radiusMinRef) > tolerance) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "radiusMin (cm): " << radiusMin << endl;
        return 4;
    }
    if (TMath::Abs(radiusMax - radiusMaxRef) > tolerance) {
        cout << "The average radius of the distribution is wrong!" << endl;
        cout << "radiusMax (cm): " << radiusMax << endl;
        return 5;
    }

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() < 1300 || run.GetEntries() > 1450) {
        cout << "The number of entries is wrong!" << endl;
        cout << "Number of entries : " << run.GetEntries() << endl;
        return 5;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";

    return 0;
}
