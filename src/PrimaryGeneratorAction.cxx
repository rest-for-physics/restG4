
#include "PrimaryGeneratorAction.h"

#include <TF1.h>
#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PrimaryGeneratorInfo.h>

#include <G4Event.hh>
#include <G4Geantino.hh>
#include <G4IonTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>

#include "SimulationManager.h"

using namespace std;
using namespace TRestGeant4PrimaryGeneratorTypes;

PrimaryGeneratorAction::PrimaryGeneratorAction(SimulationManager* simulationManager)
    : G4VUserPrimaryGeneratorAction(), fSimulationManager(simulationManager) {
    fGeneratorSpatialDensityFunction = nullptr;

    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    TRestGeant4ParticleSource* source = restG4Metadata->GetParticleSource(0);

    const string angularDistTypeName = source->GetAngularDistributionType().Data();
    const auto angularDistTypeEnum = StringToAngularDistributionTypes(angularDistTypeName);

    const string energyDistTypeName = source->GetEnergyDistributionType().Data();
    const auto energyDistTypeEnum = StringToEnergyDistributionTypes(energyDistTypeName);

    fRandom = new TRandom(restG4Metadata->GetSeed() + TRandom(G4Threading::G4GetThreadId()).Integer(1E9));

    if (energyDistTypeEnum == EnergyDistributionTypes::TH1D) {
        Double_t minEnergy = source->GetEnergyDistributionRangeMin();
        Double_t maxEnergy = source->GetEnergyDistributionRangeMax();
        SetEnergyDistributionHistogram(fSimulationManager->GetPrimaryEnergyDistribution(), minEnergy,
                                       maxEnergy);
    } else if (energyDistTypeEnum == EnergyDistributionTypes::FORMULA) {
        fEnergyDistributionFunction = (TF1*)source->GetEnergyDistributionFunction()->Clone();
        auto newRangeXMin = fEnergyDistributionFunction->GetXmin();
        if (source->GetEnergyDistributionRangeMin() > fEnergyDistributionFunction->GetXmin()) {
            newRangeXMin = source->GetEnergyDistributionRangeMin();
        }
        auto newRangeXMax = fEnergyDistributionFunction->GetXmax();
        if (source->GetEnergyDistributionRangeMax() < fEnergyDistributionFunction->GetXmax()) {
            newRangeXMax = source->GetEnergyDistributionRangeMax();
        }
        if (newRangeXMin == newRangeXMax || newRangeXMin > newRangeXMax) {
            cout << "PrimaryGeneratorAction - ERROR: energy distribution range is invalid" << endl;
            exit(1);
        }
        fEnergyDistributionFunction->SetRange(newRangeXMin, newRangeXMax);
    }

    if (angularDistTypeEnum == AngularDistributionTypes::TH1D) {
        SetAngularDistributionHistogram(fSimulationManager->GetPrimaryAngularDistribution());
    } else if (angularDistTypeEnum == AngularDistributionTypes::FORMULA) {
        fAngularDistributionFunction = (TF1*)source->GetAngularDistributionFunction()->Clone();
    }
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

void PrimaryGeneratorAction::SetEnergyDistributionHistogram(const TH1D* h, double eMin, double eMax) {
    auto xLabel = (TString)h->GetXaxis()->GetTitle();

    if (xLabel.Contains("MeV")) {
        energyFactor = 1.e3;
    } else if (xLabel.Contains("GeV")) {
        energyFactor = 1.e6;
    } else {
        energyFactor = 1.;
    }

    fEnergyDistributionHistogram = h;
    fSpectrumIntegral = fEnergyDistributionHistogram->Integral();

    startEnergyBin = 1;
    endEnergyBin = fEnergyDistributionHistogram->GetNbinsX();

    if (eMin > 0) {
        for (int i = startEnergyBin; i <= endEnergyBin; i++) {
            if (fEnergyDistributionHistogram->GetBinCenter(i) > eMin) {
                startEnergyBin = i;
                break;
            }
        }
    }

    if (eMax > 0) {
        for (int i = startEnergyBin; i <= endEnergyBin; i++) {
            if (fEnergyDistributionHistogram->GetBinCenter(i) > eMax) {
                endEnergyBin = i;
                break;
            }
        }
    }

    fSpectrumIntegral = fEnergyDistributionHistogram->Integral(startEnergyBin, endEnergyBin);
}

void PrimaryGeneratorAction::SetGeneratorSpatialDensity(TString str) {
    auto expression = (string)str;
    delete fGeneratorSpatialDensityFunction;
    if (expression.find_first_of("xyz") == -1) {
        fGeneratorSpatialDensityFunction = nullptr;
        return;
    }
    fGeneratorSpatialDensityFunction = new TF3("GeneratorDistFunc", str);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    std::lock_guard<std::mutex> lock(fMutex);  // TODO: remove this lock after fixing problems

    auto simulationManager = fSimulationManager;
    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Primary generation" << endl;
    }
    // We have to initialize here and not in start of the event because
    // GeneratePrimaries is called first, and we want to store event origin and
    // position inside
    // we should have already written the information from previous event to disk
    // (in endOfEventAction)

    for (int i = 0; i < restG4Metadata->GetNumberOfSources(); i++) {
        restG4Metadata->GetParticleSource(i)->Update();
    }

    Int_t nParticles = restG4Metadata->GetNumberOfSources();

    // Set the particle(s)' position, multiple particles generated from multiple
    // sources shall always have a same origin
    SetParticlePosition();

    for (int i = 0; i < restG4Metadata->GetNumberOfSources(); i++) {
        vector<TRestGeant4Particle> particles = restG4Metadata->GetParticleSource(i)->GetParticles();
        for (const auto& p : particles) {
            // ParticleDefinition should be always declared first (after position).
            SetParticleDefinition(i, p);
            // Particle Direction must be always set before energy
            SetParticleEnergy(i, p);
            SetParticleDirection(i, p);
            fParticleGun.GeneratePrimaryVertex(event);
        }
    }
}

G4ParticleDefinition* PrimaryGeneratorAction::SetParticleDefinition(Int_t particleSourceIndex,
                                                                    const TRestGeant4Particle& particle) {
    auto simulationManager = fSimulationManager;
    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();

    auto particleName = (string)particle.GetParticleName();

    Double_t excitedEnergy = (double)particle.GetExcitationLevel();  // in keV

    Int_t charge = particle.GetParticleCharge();

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Particle name: " << particleName << endl;
        cout << "DEBUG: Particle excited energy: " << excitedEnergy << " keV" << endl;
    }

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    if (!fParticle) {
        fParticle = particleTable->FindParticle(particleName);

        // There might be a better way to do this
        for (int Z = 1; Z <= 110; Z++) {
            for (int A = 2 * Z - 1; A <= 3 * Z; A++) {
                if (particleName == G4IonTable::GetIonTable()->GetIonName(Z, A)) {
                    // excited energy is in rest units keV, when input to geant4, we shall convert to MeV
                    fParticle = G4IonTable::GetIonTable()->GetIon(Z, A, excitedEnergy / 1000);
                    particleName = G4IonTable::GetIonTable()->GetIonName(Z, A, excitedEnergy / 1000);
                    fParticleGun.SetParticleCharge(charge);
                }
            }
        }

        if (!fParticle) {
            cout << "Particle definition : " << particleName << " not found!" << endl;
            exit(1);
        }
    }

    fParticleGun.SetParticleDefinition(fParticle);

    return fParticle;
}

void PrimaryGeneratorAction::SetParticleDirection(Int_t particleSourceIndex,
                                                  const TRestGeant4Particle& particle) {
    auto simulationManager = fSimulationManager;
    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();
    TRestGeant4ParticleSource* source = restG4Metadata->GetParticleSource(0);

    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    G4ThreeVector direction = {source->GetDirection().X(), source->GetDirection().Y(),
                               source->GetDirection().Z()};

    const string angularDistTypeName = source->GetAngularDistributionType().Data();
    const auto angularDistTypeEnum = StringToAngularDistributionTypes(angularDistTypeName);

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Angular distribution: " << angularDistTypeName << endl;
    }

    // generator type
    const string& spatialGeneratorTypeName = primaryGeneratorInfo.GetSpatialGeneratorType().Data();
    const auto spatialGeneratorTypeEnum = StringToSpatialGeneratorTypes(spatialGeneratorTypeName);

    const string& spatialGeneratorShapeName = primaryGeneratorInfo.GetSpatialGeneratorShape().Data();
    const auto spatialGeneratorShapeEnum = StringToSpatialGeneratorShapes(spatialGeneratorShapeName);

    if (angularDistTypeEnum == AngularDistributionTypes::ISOTROPIC) {
        direction = GetIsotropicVector();
    } else if (angularDistTypeEnum == AngularDistributionTypes::TH1D) {
        Double_t angle = 0;
        Double_t value = G4UniformRand() * fAngularDistributionHistogram->Integral();
        Double_t sum = 0;
        // deltaAngle is the constant x distance between bins
        Double_t deltaAngle =
            fAngularDistributionHistogram->GetBinCenter(2) - fAngularDistributionHistogram->GetBinCenter(1);
        // we sample the CDF (uniform between 0 and the distribution integral which should be equal to 1)
        // the inverse of CDF of the uniformly sampled value will follow a distribution given by the PDF, we
        // compute this inverse
        for (int bin = 1; bin <= fAngularDistributionHistogram->GetNbinsX(); bin++) {
            sum += fAngularDistributionHistogram->GetBinContent(bin);
            if (sum >= value) {
                angle =
                    fAngularDistributionHistogram->GetBinCenter(bin) + deltaAngle * (0.5 - G4UniformRand());
                break;
            }
        }

        G4ThreeVector referenceOrigin = direction;

        // We generate the distribution angle (theta) using a rotation around the orthogonal vector
        direction.rotate(angle, direction.orthogonal());

        // We rotate a full-2PI random angle along the original direction to generate a cone
        direction.rotate(G4UniformRand() * 2 * M_PI, referenceOrigin);

    } else if (angularDistTypeEnum == AngularDistributionTypes::FORMULA) {
        G4ThreeVector referenceOrigin = direction;

        // We generate the distribution angle (theta) using a rotation around the orthogonal vector
        direction.rotate(fAngularDistributionFunction->GetRandom(fRandom), direction.orthogonal());

        // We rotate a full-2PI random angle along the original direction to generate a cone
        direction.rotate(G4UniformRand() * 2 * M_PI, referenceOrigin);

    } else if (angularDistTypeEnum == AngularDistributionTypes::FLUX) {
        const TVector3& v = particle.GetMomentumDirection().Unit();
        direction.set(v.X(), v.Y(), v.Z());

    } else if (angularDistTypeEnum == AngularDistributionTypes::BACK_TO_BACK) {
        // This should never crash. In TRestG4Metadata we have defined that if the
        // first source is back to back we set it to isotropic
        // TVector3 v = restG4Event->GetPrimaryEventDirection(particleSourceIndex - 1);
        // v = v.Unit();

        // direction.set(-v.X(), -v.Y(), -v.Z());
        exit(1);
    } else {
        G4cout << "WARNING: Generator angular distribution was not recognized. Particle direction set to ("
               << direction.x() << ", " << direction.y() << ", " << direction.z() << ")" << G4endl;
    }

    fParticleGun.SetParticleMomentumDirection(direction);
}

void PrimaryGeneratorAction::SetParticleEnergy(Int_t particleSourceIndex,
                                               const TRestGeant4Particle& particle) {
    auto simulationManager = fSimulationManager;

    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();
    TRestGeant4ParticleSource* source = restG4Metadata->GetParticleSource(0);
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const string angularDistTypeName =
        restG4Metadata->GetParticleSource(particleSourceIndex)->GetAngularDistributionType().Data();
    const auto angularDistTypeEnum = StringToAngularDistributionTypes(angularDistTypeName);

    const string energyDistTypeName =
        restG4Metadata->GetParticleSource(particleSourceIndex)->GetEnergyDistributionType().Data();
    const auto energyDistTypeEnum = StringToEnergyDistributionTypes(energyDistTypeName);

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Energy distribution: " << energyDistTypeName << endl;
    }

    Double_t energy = 1 * keV;

    if (energyDistTypeEnum == EnergyDistributionTypes::MONO) {
        energy = particle.GetEnergy() * keV;
    } else if (energyDistTypeEnum == EnergyDistributionTypes::FLAT) {
        TVector2 enRange =
            restG4Metadata->GetParticleSource(particleSourceIndex)->GetEnergyDistributionRange();
        energy = ((enRange.Y() - enRange.X()) * G4UniformRand() + enRange.X()) * keV;
    } else if (energyDistTypeEnum == EnergyDistributionTypes::LOG) {
        TVector2 enRange =
            restG4Metadata->GetParticleSource(particleSourceIndex)->GetEnergyDistributionRange();
        auto max_energy = enRange.Y() * keV;
        auto min_energy = enRange.X() * keV;
        energy = exp((log(max_energy) - log(min_energy)) * G4UniformRand() + log(min_energy));

    } else if (energyDistTypeEnum == EnergyDistributionTypes::TH1D) {
        Double_t value = G4UniformRand() * fSpectrumIntegral;
        Double_t sum = 0;
        Double_t deltaEnergy =
            fEnergyDistributionHistogram->GetBinCenter(2) - fEnergyDistributionHistogram->GetBinCenter(1);
        for (int bin = startEnergyBin; bin <= endEnergyBin; bin++) {
            sum += fEnergyDistributionHistogram->GetBinContent(bin);

            if (sum > value) {
                energy = energyFactor *
                         (Double_t)(fEnergyDistributionHistogram->GetBinCenter(bin) +
                                    deltaEnergy * (0.5 - G4UniformRand())) *
                         keV;
                break;
            }
        }
    } else if (energyDistTypeEnum == EnergyDistributionTypes::FORMULA) {
        energy = fEnergyDistributionFunction->GetRandom(fRandom) * keV;
    } else {
        G4cout << "WARNING! Energy distribution type was not recognized. Setting "
                  "energy to 1keV"
               << G4endl;
        energy = 1 * keV;
    }

    if (particleSourceIndex > 0 &&
        angularDistTypeEnum == TRestGeant4PrimaryGeneratorTypes::AngularDistributionTypes::BACK_TO_BACK) {
        energy = lastEnergy;
    }

    if (particleSourceIndex == 0) {
        lastEnergy = energy;
    }
    fParticleGun.SetParticleEnergy(energy);

    if (restG4Metadata->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Particle energy: " << energy / keV << " keV" << endl;
    }
}

void PrimaryGeneratorAction::SetParticlePosition() {
    auto simulationManager = fSimulationManager;

    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const string& spatialGeneratorTypeName = primaryGeneratorInfo.GetSpatialGeneratorType().Data();
    const auto spatialGeneratorTypeEnum = StringToSpatialGeneratorTypes(spatialGeneratorTypeName);

    const string& spatialGeneratorShapeName = primaryGeneratorInfo.GetSpatialGeneratorShape().Data();
    const auto spatialGeneratorShapeEnum = StringToSpatialGeneratorShapes(spatialGeneratorShapeName);

    double x = 0, y = 0, z = 0;

    while (true) {
        if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::POINT) {
            GenPositionOnPoint(x, y, z);
        } else if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::SURFACE) {
            if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::GDML) {
                GenPositionOnGDMLSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::BOX) {
                GenPositionOnBoxSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::CYLINDER) {
                GenPositionOnCylinderSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::SPHERE) {
                GenPositionOnSphereSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::CIRCLE) {
                GenPositionOnDisk(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::WALL) {
                GenPositionOnWall(x, y, z);
            }
        } else if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::VOLUME) {
            if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::GDML) {
                GenPositionOnGDMLVolume(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::BOX) {
                GenPositionOnBoxVolume(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::CYLINDER) {
                GenPositionOnCylinderVolume(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::SPHERE) {
                GenPositionOnSphereVolume(x, y, z);
            }
        } else {
            G4cout << "WARNING! Generator type \"" << spatialGeneratorTypeName
                   << "\" was not recognized. Launching particle "
                      "from origin (0,0,0)"
                   << G4endl;
        }

        // use the density function. If the density is small, then val2 is small, we are more
        // likely to regenerate the particle position
        if (fGeneratorSpatialDensityFunction) {
            double val1 = G4UniformRand();
            double val2 = fGeneratorSpatialDensityFunction->Eval(x, y, z);
            if (val2 > 1) {
                cout << "error! Generator density function > 1 at position (" << x << ", " << y << ", " << z
                     << "), check your definition!" << endl;
                exit(1);
            }
            if (val1 > val2) {
                continue;
            }
        }
        break;
    }

    fParticleGun.SetParticlePosition(G4ThreeVector(x, y, z));
}

G4ThreeVector PrimaryGeneratorAction::GetIsotropicVector() const {
    G4double a, b, c;
    G4double n;
    do {
        a = (G4UniformRand() - 0.5) / 0.5;
        b = (G4UniformRand() - 0.5) / 0.5;
        c = (G4UniformRand() - 0.5) / 0.5;
        n = a * a + b * b + c * c;
    } while (n > 1 || n == 0.0);

    n = sqrt(n);
    a /= n;
    b /= n;
    c /= n;
    return {a, b, c};
}

Double_t PrimaryGeneratorAction::GetAngle(G4ThreeVector x, G4ThreeVector y) {
    Double_t angle = y.angle(x);

    return angle;
}

Double_t PrimaryGeneratorAction::GetCosineLowRandomThetaAngle() {
    // We obtain an angle with a cos(theta)*sin(theta) distribution
    Double_t value = G4UniformRand();
    double dTheta = 0.01;
    for (double theta = 0; theta < M_PI / 2; theta = theta + dTheta) {
        // sin(theta)^2 is the integral of sin(theta)*cos(theta)
        if (sin(theta) * sin(theta) >= value) {
            if (theta > 0)
                return theta - dTheta * (0.5 - G4UniformRand());
            else
                return theta;
        }
    }
    return M_PI / 2;
}

void PrimaryGeneratorAction::GenPositionOnGDMLVolume(double& x, double& y, double& z) {
    auto detector = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    double xMin = detector->GetBoundBoxXMin();
    double xMax = detector->GetBoundBoxXMax();
    double yMin = detector->GetBoundBoxYMin();
    double yMax = detector->GetBoundBoxYMax();
    double zMin = detector->GetBoundBoxZMin();
    double zMax = detector->GetBoundBoxZMax();

    do {
        x = xMin + (xMax - xMin) * G4UniformRand();
        y = yMin + (yMax - yMin) * G4UniformRand();
        z = zMin + (zMax - zMin) * G4UniformRand();
    } while (detector->GetGeneratorSolid()->Inside(G4ThreeVector(x, y, z)) != kInside);

    x = x + detector->GetGeneratorTranslation().x();
    y = y + detector->GetGeneratorTranslation().y();
    z = z + detector->GetGeneratorTranslation().z();
}

void PrimaryGeneratorAction::GenPositionOnGDMLSurface(double& x, double& y, double& z) {
    // TODO there is a problem, probably with G4 function GetPointOnSurface
    // It produces a point on the surface but it is not uniformly distributed
    // May be it is just an OPENGL drawing issue?
    auto detector = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    G4ThreeVector position = detector->GetGeneratorSolid()->GetPointOnSurface();

    x = position.x();
    y = position.y();
    z = position.z();

    x = x + detector->GetGeneratorTranslation().x();
    y = y + detector->GetGeneratorTranslation().y();
    z = z + detector->GetGeneratorTranslation().z();
}

void PrimaryGeneratorAction::GenPositionOnBoxVolume(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    Double_t sizeX = primaryGeneratorInfo.GetSpatialGeneratorSize().X();
    Double_t sizeY = primaryGeneratorInfo.GetSpatialGeneratorSize().Y();
    Double_t sizeZ = primaryGeneratorInfo.GetSpatialGeneratorSize().Z();

    x = sizeX * (G4UniformRand() - 0.5);
    y = sizeY * (G4UniformRand() - 0.5);
    z = sizeZ * (G4UniformRand() - 0.5);

    G4ThreeVector position = G4ThreeVector(x, y, z);

    const TVector3& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const TVector3& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnBoxSurface(double& x, double& y, double& z) {
    cout << __PRETTY_FUNCTION__ << ": not implemented!" << endl;
    exit(1);
}

void PrimaryGeneratorAction::GenPositionOnSphereVolume(double& x, double& y, double& z) {
    cout << __PRETTY_FUNCTION__ << ": not implemented!" << endl;
    exit(1);
}

void PrimaryGeneratorAction::GenPositionOnSphereSurface(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    G4ThreeVector position = GetIsotropicVector();

    const Double_t radius = primaryGeneratorInfo.GetSpatialGeneratorSize().X();

    const TVector3& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = radius * position.x() + center.X();
    y = radius * position.y() + center.Y();
    z = radius * position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnCylinderVolume(double& x, double& y, double& z) {
    cout << __PRETTY_FUNCTION__ << ": not implemented!" << endl;
    exit(1);
}

void PrimaryGeneratorAction::GenPositionOnCylinderSurface(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    Double_t angle = 2 * M_PI * G4UniformRand();

    const Double_t radius = primaryGeneratorInfo.GetSpatialGeneratorSize().X();
    const Double_t length = primaryGeneratorInfo.GetSpatialGeneratorSize().Y();

    x = radius * cos(angle);
    y = radius * sin(angle);
    z = length * (G4UniformRand() - 0.5);

    G4ThreeVector position = G4ThreeVector(x, y, z);

    const TVector3& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const TVector3& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnPoint(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const TVector3& position = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.X();
    y = position.Y();
    z = position.Z();
}

void PrimaryGeneratorAction::GenPositionOnWall(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const Double_t sizeX = primaryGeneratorInfo.GetSpatialGeneratorSize().X();
    const Double_t sizeY = primaryGeneratorInfo.GetSpatialGeneratorSize().Y();

    x = sizeX * (G4UniformRand() - 0.5);
    y = sizeY * (G4UniformRand() - 0.5);

    G4ThreeVector position = G4ThreeVector(x, y, 0);

    const TVector3& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const TVector3& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnDisk(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const Double_t radius = primaryGeneratorInfo.GetSpatialGeneratorSize().X();

    do {
        x = 2 * radius * (G4UniformRand() - 0.5);
        y = 2 * radius * (G4UniformRand() - 0.5);
    } while (x * x + y * y > radius * radius);

    G4ThreeVector position = G4ThreeVector(x, y, 0);

    const TVector3& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const TVector3& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}
