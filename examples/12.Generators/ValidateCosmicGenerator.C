#include <TRestGeant4Event.h>

Int_t ValidateCosmicGenerator(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);

    cout << "Run entries: " << run.GetEntries() << endl;
    if (run.GetEntries() != 1E6) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 2;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    constexpr double tolerance = 0.1;

    TVector3 sourceDirection = metadata->GetParticleSource()->GetDirection();
    if (TVector3(0.1, -0.5, 1.25).Angle(sourceDirection) > 0.01) {
        cout << "Bad source direction: (" << sourceDirection.X() << ", " << sourceDirection.Y() << ", "
             << sourceDirection.Z() << endl;
        return 3;
    }

    double thetaAverage = 0, thetaMin = TMath::Infinity(), thetaMax = 0;
    constexpr double thetaAverageRef = 37.8595, thetaMinRef = 0.0110538, thetaMaxRef = 89.9982;

    double energyPrimaryAverage = 0, energyPrimaryMin = TMath::Infinity(), energyPrimaryMax = 0;
    constexpr double energyPrimaryAverageRef = 7.89997e+06, energyPrimaryMinRef = 200008,
                     energyPrimaryMaxRef = 4.63073e+09;

    double primaryRadiusRef = 173.20508;

    TH1D thetaHist("thetaHist", "Theta angle from source direction", 100, 0, TMath::Pi() * TMath::RadToDeg());
    TH1D energyHist("energyHist", "Primary muon energy", 200, 0, 1E9);
    for (int i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);
        TVector3 direction = event->GetPrimaryEventDirection();
        TVector3 position = event->GetPrimaryEventOrigin();
        if (TMath::Abs(position.Mag() - primaryRadiusRef) / primaryRadiusRef > tolerance) {
            cout << "Bad primary position: (" << position.X() << ", " << position.Y() << ", " << position.Z()
                 << ") with radius " << position.Mag() << endl;
            return 4;
        }
        const auto theta = sourceDirection.Angle(direction) * TMath::RadToDeg();
        thetaHist.Fill(theta);

        Double_t energy = event->GetPrimaryEventEnergy();
        energyHist.Fill(energy);

        thetaAverage += theta / run.GetEntries();
        if (theta < thetaMin) {
            thetaMin = theta;
        }
        if (theta > thetaMax) {
            thetaMax = theta;
        }

        energyPrimaryAverage += energy / run.GetEntries();
        if (energy < energyPrimaryMin) {
            energyPrimaryMin = energy;
        }
        if (energy > energyPrimaryMax) {
            energyPrimaryMax = energy;
        }
    }
    thetaHist.Scale(1. / thetaHist.Integral(), "width");
    energyHist.Scale(1. / energyHist.Integral(), "width");

    cout << "Average theta (deg): " << thetaAverage << endl;
    cout << "Minimum theta (deg): " << thetaMin << endl;
    cout << "Maximum theta (deg): " << thetaMax << endl;

    cout << "Average energy (keV): " << energyPrimaryAverage << endl;
    cout << "Minimum energy (keV): " << energyPrimaryMin << endl;
    cout << "Maximum energy (keV): " << energyPrimaryMax << endl;

    if (TMath::Abs(thetaAverage - thetaAverageRef) / thetaAverageRef > tolerance) {
        cout << "The average theta of the distribution is wrong!" << endl;
        cout << "thetaAverage (deg): " << thetaAverage << endl;
        return 5;
    }
    if (TMath::Abs(thetaMin - thetaMinRef) / (thetaMinRef + 1) > tolerance) {
        cout << "The minimum theta of the distribution is wrong!" << endl;
        cout << "thetaMin (deg): " << thetaMin << endl;
        return 6;
    }
    if (TMath::Abs(thetaMax - thetaMaxRef) / thetaMaxRef > tolerance) {
        cout << "The maximum theta of the distribution is wrong!" << endl;
        cout << "thetaMax (deg): " << thetaMax << endl;
        return 7;
    }

    if (TMath::Abs(energyPrimaryAverage - energyPrimaryAverageRef) / energyPrimaryAverageRef > tolerance) {
        cout << "The average energy of the distribution is wrong!" << endl;
        cout << "energyPrimaryAverage (keV): " << energyPrimaryAverage << endl;
        return 8;
    }
    if (TMath::Abs(energyPrimaryMin - energyPrimaryMinRef) / energyPrimaryMinRef > tolerance) {
        cout << "The minimum energy of the distribution is wrong!" << endl;
        cout << "energyPrimaryMin (keV): " << energyPrimaryMin << endl;
        return 9;
    }
    if (TMath::Abs(energyPrimaryMax - energyPrimaryMaxRef) / energyPrimaryMaxRef > tolerance) {
        cout << "The maximum energy of the distribution is wrong!" << endl;
        cout << "energyPrimaryMax (keV): " << energyPrimaryMax << endl;
        return 10;
    }

    if (metadata->GetNumberOfEvents() != 1E6) {
        cout << "wrong number of events: " << metadata->GetNumberOfEvents() << endl;
        return 11;
    }

    const auto surfaceTermRef = 942.4778;
    if (abs(metadata->GetGeant4PrimaryGeneratorInfo().GetSpatialGeneratorCosmicSurfaceTermCm2() -
            surfaceTermRef) /
            surfaceTermRef >
        tolerance) {
        cout << "wrong cosmic surface term: "
             << metadata->GetGeant4PrimaryGeneratorInfo().GetSpatialGeneratorCosmicSurfaceTermCm2() << endl;
        return 12;
    }

    const auto cosmicFluxRef = 0.108215;
    if (abs(metadata->GetCosmicFluxInCountsPerCm2PerSecond() - cosmicFluxRef) / cosmicFluxRef > tolerance) {
        cout << "wrong cosmic flux: " << metadata->GetCosmicFluxInCountsPerCm2PerSecond() << endl;
        return 13;
    }

    const auto simulationTimeRef = 9804.87;
    if (abs(metadata->GetEquivalentSimulatedTime() - simulationTimeRef) / simulationTimeRef > tolerance) {
        cout << "wrong equivalent simulation time: " << metadata->GetEquivalentSimulatedTime() << endl;
        return 14;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
