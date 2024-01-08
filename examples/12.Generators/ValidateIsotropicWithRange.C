#include <TRestGeant4Event.h>

Int_t ValidateEnergyAndAngularRange(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() != 1000000) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 2;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    TVector3 sourceDirection = metadata->GetParticleSource()->GetDirection();
    cout << "original source direction: (" << sourceDirection.X() << ", " << sourceDirection.Y() << ", "
         << sourceDirection.Z() << ")" << endl;

    constexpr double tolerance = 0.1;

    double thetaAverage = 0, thetaMin = TMath::Infinity(), thetaMax = 0;
    constexpr double thetaAverageRef = 23.1283, thetaMinRef = 0.0, thetaMaxRef = 55.0;

    double energyPrimaryAverage = 0, energyPrimaryMin = TMath::Infinity(), energyPrimaryMax = 0;
    constexpr double energyPrimaryAverageRef = 2915.84, energyPrimaryMinRef = 100,
                     energyPrimaryMaxRef = 15000.0;

    for (int i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);
        TVector3 direction = event->GetPrimaryEventDirection();
        const auto theta = sourceDirection.Angle(direction) * TMath::RadToDeg();

        thetaAverage += theta / run.GetEntries();
        if (theta < thetaMin) {
            thetaMin = theta;
        }
        if (theta > thetaMax) {
            thetaMax = theta;
        }

        Double_t energy = event->GetPrimaryEventEnergy();
        energyPrimaryAverage += energy / run.GetEntries();
        if (energy < energyPrimaryMin) {
            energyPrimaryMin = energy;
        }
        if (energy > energyPrimaryMax) {
            energyPrimaryMax = energy;
        }
    }

    cout << "Average theta (deg): " << thetaAverage << endl;
    cout << "Minimum theta (deg): " << thetaMin << endl;
    cout << "Maximum theta (deg): " << thetaMax << endl;

    cout << "Average energy (keV): " << energyPrimaryAverage << endl;
    cout << "Minimum energy (keV): " << energyPrimaryMin << endl;
    cout << "Maximum energy (keV): " << energyPrimaryMax << endl;

    if (TMath::Abs(thetaAverage - thetaAverageRef) / thetaAverageRef > tolerance) {
        cout << "The average theta of the distribution is wrong!" << endl;
        cout << "thetaAverage (deg): " << thetaAverage << endl;
        return 3;
    }
    if (TMath::Abs(thetaMin - thetaMinRef) / (thetaMinRef + 1) > tolerance) {
        cout << "The minimum theta of the distribution is wrong!" << endl;
        cout << "thetaMin (deg): " << thetaMin << endl;
        return 4;
    }
    if (TMath::Abs(thetaMax - thetaMaxRef) / thetaMaxRef > tolerance) {
        cout << "The maximum theta of the distribution is wrong!" << endl;
        cout << "thetaMax (deg): " << thetaMax << endl;
        return 5;
    }

    if (TMath::Abs(energyPrimaryAverage - energyPrimaryAverageRef) / energyPrimaryAverageRef > tolerance) {
        cout << "The average energy of the distribution is wrong!" << endl;
        cout << "energyPrimaryAverage (keV): " << energyPrimaryAverage << endl;
        return 6;
    }
    if (TMath::Abs(energyPrimaryMin - energyPrimaryMinRef) / energyPrimaryMinRef > tolerance) {
        cout << "The minimum energy of the distribution is wrong!" << endl;
        cout << "energyPrimaryMin (keV): " << energyPrimaryMin << endl;
        return 7;
    }
    if (TMath::Abs(energyPrimaryMax - energyPrimaryMaxRef) / energyPrimaryMaxRef > tolerance) {
        cout << "The maximum energy of the distribution is wrong!" << endl;
        cout << "energyPrimaryMax (keV): " << energyPrimaryMax << endl;
        return 8;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
