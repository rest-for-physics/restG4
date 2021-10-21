
#include "PrimaryGeneratorAction.h"

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <spdlog/spdlog.h>

#include <G4Event.hh>
#include <G4Geantino.hh>
#include <G4IonTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>

#include "GlobalManager.h"

Int_t face = 0;

double GeneratorRndm() { return G4UniformRand(); }

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(make_unique<G4ParticleGun>()),
      fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()) {
    spdlog::info("PrimaryGeneratorAction::PrimaryGeneratorAction");

    TH1D* primaryEnergyDistribution = GlobalManager::Instance()->GetPrimaryEnergyDistribution();
    if (primaryEnergyDistribution) {
        SetEnergySpectrum(*primaryEnergyDistribution, 0, 0);
        spdlog::info(
            "PrimaryGeneratorAction:PrimaryGeneratorAction - Setting primary energy distribution to {}",
            fPrimaryEnergySpectrum.GetName());
    }

    auto primaryAngularDistribution = GlobalManager::Instance()->GetPrimaryAngularDistribution();
    if (primaryAngularDistribution) {
        SetAngularDistribution(*primaryAngularDistribution);
        spdlog::info(
            "PrimaryGeneratorAction:PrimaryGeneratorAction - Setting primary angular distribution to {}",
            fPrimaryAngularSpectrum.GetName());
    }

    fDetector = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    G4int n_particle = 1;
    fParticleGun->SetNumberOfParticles(n_particle);

    fGeneratorSpatialDensityFunction = nullptr;

    for (int i = 0; i < fRestGeant4Metadata->GetNumberOfSources(); i++) {
        fRestGeant4Metadata->GetParticleSource(i)->SetRndmMethod(GeneratorRndm);
    }
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

void PrimaryGeneratorAction::SetEnergySpectrum(const TH1D& h, double eMin, double eMax) {
    TString xLabel = h.GetXaxis()->GetTitle();

    if (xLabel.Contains("MeV")) {
        energyFactor = 1.e3;
    } else if (xLabel.Contains("GeV")) {
        energyFactor = 1.e6;
    } else {
        energyFactor = 1.;
    }

    fPrimaryEnergySpectrum = h;
    fSpectrumIntegral = fPrimaryEnergySpectrum.Integral();

    startEnergyBin = 1;
    endEnergyBin = fPrimaryEnergySpectrum.GetNbinsX();

    if (eMin > 0) {
        for (int i = startEnergyBin; i <= endEnergyBin; i++) {
            if (fPrimaryEnergySpectrum.GetBinCenter(i) > eMin) {
                startEnergyBin = i;
                break;
            }
        }
    }

    if (eMax > 0) {
        for (int i = startEnergyBin; i <= endEnergyBin; i++) {
            if (fPrimaryEnergySpectrum.GetBinCenter(i) > eMax) {
                endEnergyBin = i;
                break;
            }
        }
    }

    fSpectrumIntegral = fPrimaryEnergySpectrum.Integral(startEnergyBin, endEnergyBin);
}

void PrimaryGeneratorAction::SetAngularDistribution(const TH1D& h) { fPrimaryAngularSpectrum = h; }

void PrimaryGeneratorAction::SetGeneratorSpatialDensity(TString str) {
    string expression = (string)str;
    delete fGeneratorSpatialDensityFunction;
    if (expression.find_first_of("xyz") == -1) {
        fGeneratorSpatialDensityFunction = nullptr;
        return;
    }
    fGeneratorSpatialDensityFunction = new TF3("GeneratorDistFunc", str);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    // return;
    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Debug) {
        G4cout << "DEBUG: Primary generation" << endl;
    }

    for (int i = 0; i < fRestGeant4Metadata->GetNumberOfSources(); i++) {
        fRestGeant4Metadata->GetParticleSource(i)->Update();
    }

    Int_t nParticles = fRestGeant4Metadata->GetNumberOfSources();

    // Set the particle(s)' position, multiple particles generated from multiple
    // sources shall always have a same origin
    SetParticlePosition();

    for (int i = 0; i < fRestGeant4Metadata->GetNumberOfSources(); i++) {
        for (const auto& particle : fRestGeant4Metadata->GetParticleSource(i)->GetParticles()) {
            // ParticleDefinition should be always declared first (after position).
            SetParticleDefinition(i, particle);
            SetParticleEnergy(i, particle);
            SetParticleDirection(i, particle);

            fParticleGun->GeneratePrimaryVertex(event);
        }
    }

    auto primary = event->GetPrimaryVertex()->GetPrimary();
    const auto primaryEnergy = primary->GetKineticEnergy();
    const auto& primaryPosition = event->GetPrimaryVertex()->GetPosition();
    const auto& primaryDirection = primary->GetMomentumDirection();

    spdlog::debug(
        "PrimaryGeneratorAction::GeneratePrimaries - Particle: {} - Position (mm): [{:03.2f}, {:03.2f}, "
        "{:03.2f}] - Energy: {:03.2f} keV"
        " - Direction: [{:03.4f}, {:03.4f}, {:03.4f}] ",                               //
        primary->GetParticleDefinition()->GetParticleName(),                           //
        primaryPosition.x() / mm, primaryPosition.y() / mm, primaryPosition.z() / mm,  //
        primaryEnergy / keV,                                                           //
        primaryDirection.x(), primaryDirection.y(), primaryDirection.z()               //
    );
}

G4ParticleDefinition* PrimaryGeneratorAction::SetParticleDefinition(Int_t n, TRestGeant4Particle p) {
    string particle_name = (string)p.GetParticleName();

    auto excited_energy = (double)p.GetExcitationLevel();  // in keV

    Int_t charge = p.GetParticleCharge();

    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Debug) {
        G4cout << "DEBUG: Particle name: " << particle_name << endl;
        // G4cout << "DEBUG: Particle charge: " << charge << endl;
        G4cout << "DEBUG: Particle excited energy: " << excited_energy << " keV" << endl;
    }

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    if (!fParticle) {
        fParticle = particleTable->FindParticle(particle_name);

        // There might be a better way to do this
        for (int Z = 1; Z <= 110; Z++)
            for (int A = 2 * Z - 1; A <= 3 * Z; A++) {
                if (particle_name == G4IonTable::GetIonTable()->GetIonName(Z, A)) {
                    // excited energy is in rest units keV, when input to geant4, we shall convert to MeV
                    fParticle = G4IonTable::GetIonTable()->GetIon(Z, A, excited_energy / 1000);
                    particle_name = G4IonTable::GetIonTable()->GetIonName(Z, A, excited_energy / 1000);
                    fParticleGun->SetParticleCharge(charge);
                }
            }
        // }

        if (!fParticle) {
            G4cout << "Particle definition : " << particle_name << " not found!" << endl;
            exit(1);
        }
    }

    fParticleGun->SetParticleDefinition(fParticle);

    return fParticle;
}

void PrimaryGeneratorAction::SetParticleDirection(Int_t n, TRestGeant4Particle p) {
    G4ThreeVector direction;
    // TODO: maybe reduce code redundancy by defining some functions?
    // TODO: fix bug when giving TH1D with lowercase (e.g. Th1D). string conversion is OK but integral gives
    // exception.
    string angular_dist_type_name = (string)fRestGeant4Metadata->GetParticleSource(n)->GetAngularDistType();
    angular_dist_type_name = g4_metadata_parameters::CleanString(angular_dist_type_name);
    g4_metadata_parameters::angular_dist_types angular_dist_type;

    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Debug) {
        G4cout << "DEBUG: Angular distribution: " << angular_dist_type_name << endl;
    }
    // we first check if it is a valid parameter
    if (g4_metadata_parameters::angular_dist_types_map.count(angular_dist_type_name)) {
        angular_dist_type = g4_metadata_parameters::angular_dist_types_map[angular_dist_type_name];
    } else {
        // if we get here it means the parameter is not valid, we can either assign a default value or stop
        // execution default value
        G4cout << "Invalid angular distribution (" + angular_dist_type_name + ") valid values are: ";
        for (const auto& pair : g4_metadata_parameters::angular_dist_types_map) {
            G4cout << pair.first << ", ";
        }
        G4cout << std::endl;
        throw "Invalid angular distribution";
    }
    // generator type
    string generator_type_name = (string)fRestGeant4Metadata->GetGeneratorType();
    generator_type_name = g4_metadata_parameters::CleanString(generator_type_name);
    g4_metadata_parameters::generator_types generator_type;
    if (g4_metadata_parameters::generator_types_map.count(generator_type_name)) {
        generator_type = g4_metadata_parameters::generator_types_map[generator_type_name];
    } else {
        // if we get here it means the parameter is not valid, we can either assign a default value or stop
        // execution default value
        G4cout << "Invalid generator type (" + generator_type_name + ") valid values are: ";
        for (const auto& pair : g4_metadata_parameters::generator_types_map) {
            G4cout << pair.first << ", ";
        }
        G4cout << std::endl;
        throw "Invalid generator type";
    }

    if (angular_dist_type == g4_metadata_parameters::angular_dist_types::ISOTROPIC) {
        // if (generator_type == g4_metadata_parameters::generator_types::VIRTUAL_BOX) {
        //    if (face == 0) direction.set(0, -1, 0);
        //    if (face == 1) direction.set(0, 1, 0);
        //    if (face == 2) direction.set(-1, 0, 0);
        //    if (face == 3) direction.set(1, 0, 0);
        //    if (face == 4) direction.set(0, 0, -1);
        //    if (face == 5) direction.set(0, 0, 1);

        //    Double_t theta = GetCosineLowRandomThetaAngle();
        //    // recording the primaries distribution
        //    G4ThreeVector referenceOrigin = direction;

        //    // We rotate the origin direction by the angular distribution angle
        //    G4ThreeVector orthoVector = direction.orthogonal();
        //    direction.rotate(theta, orthoVector);

        //    // We rotate a random angle along the original direction
        //    Double_t randomAngle = G4UniformRand() * 2 * M_PI;
        //    direction.rotate(randomAngle, referenceOrigin);
        //} else if (generator_type == g4_metadata_parameters::generator_types::VIRTUAL_SPHERE) {
        //    direction = -fParticleGun->GetParticlePosition().unit();

        //    Double_t theta = GetCosineLowRandomThetaAngle();

        //    G4ThreeVector referenceOrigin = direction;

        //    // We rotate the origin direction by the angular distribution angle
        //    G4ThreeVector orthoVector = direction.orthogonal();
        //    direction.rotate(theta, orthoVector);

        //    // We rotate a random angle along the original direction
        //    Double_t randomAngle = G4UniformRand() * 2 * M_PI;
        //    direction.rotate(randomAngle, referenceOrigin);

        //} else {
        direction = GetIsotropicVector();
        //}
    } else if (angular_dist_type == g4_metadata_parameters::angular_dist_types::TH1D) {
        Double_t angle = 0;
        Double_t value = G4UniformRand() * (fPrimaryAngularSpectrum.Integral());
        Double_t sum = 0;
        // deltaAngle is the constant x distance between bins
        Double_t deltaAngle =
            fPrimaryAngularSpectrum.GetBinCenter(2) - fPrimaryAngularSpectrum.GetBinCenter(1);
        // we sample the CDF (uniform between 0 and the distribution integral which should be equal to 1)
        // the inverse of CDF of the uniformly sampled value will follow a distribution given by the PDF, we
        // compute this inverse
        for (int bin = 1; bin <= fPrimaryAngularSpectrum.GetNbinsX(); bin++) {
            sum += fPrimaryAngularSpectrum.GetBinContent(bin);

            if (sum >= value) {
                angle = fPrimaryAngularSpectrum.GetBinCenter(bin) + deltaAngle * (0.5 - G4UniformRand());
                break;
            }
        }

        // A vector pointing to the origin (virtualSphere )
        direction = -fParticleGun->GetParticlePosition().unit();

        if (direction.x() == 0 && direction.y() == 0 && direction.z() == 0) {
            G4cout << "----------------------------------------------------------------"
                    "-----"
                 << endl;
            G4cout << "REST WARNNING : Particle being launched from the ORIGIN!! Wrong "
                    "momentum direction!"
                 << endl;
            G4cout << "Setting direction to (1,0,0)" << endl;
            G4cout << "REST angular distribution is just implemented for virtualBox "
                    "and virtualSphere"
                 << endl;
            G4cout << "Other spatial distributions can be set but it will launch the "
                    "event\n with a distribution direction to the origin of "
                    "coordinates"
                 << endl;
            G4cout << "----------------------------------------------------------------"
                    "-----"
                 << endl;
            direction.set(1, 0, 0);
        }

        // if (generator_type == g4_metadata_parameters::generator_types::VIRTUAL_BOX) {
        //    if (face == 0) direction.set(0, -1, 0);
        //    if (face == 1) direction.set(0, 1, 0);
        //    if (face == 2) direction.set(-1, 0, 0);
        //    if (face == 3) direction.set(1, 0, 0);
        //    if (face == 4) direction.set(0, 0, -1);
        //    if (face == 5) direction.set(0, 0, 1);
        //}

        // if (generator_type == g4_metadata_parameters::generator_types::VIRTUAL_WALL) {
        //    /*
        //    The default plane (virtualWall) is an XY plane so the default normal vector is (0,0,1).
        //    We will rotate this vector according to the generator rotation so that keeps being normal to the
        //    plane
        //    */
        //    TVector3 normal(0, 0, 1);

        //    normal.RotateX(M_PI * fRestGeant4Metadata->GetGeneratorRotation().X() / 180);
        //    normal.RotateY(M_PI * fRestGeant4Metadata->GetGeneratorRotation().Y() / 180);
        //    normal.RotateZ(M_PI * fRestGeant4Metadata->GetGeneratorRotation().Z() / 180);

        //    /*
        //    Depending on which rotation we chose for the plane the normal vector can now point outwards of
        //    the detector. We rotate so that it always looks towards the center (0,0,0)
        //    */
        //    TVector3 generator_position = fRestGeant4Metadata->GetGeneratorPosition().Unit();
        //    if (generator_position.x() * normal.x() + generator_position.y() * normal.y() +
        //            generator_position.z() * normal.z() >
        //        0) {
        //        normal = (-1) * normal;
        //    }

        //    direction.set(normal.x(), normal.y(), normal.z());
        //}

        G4ThreeVector referenceOrigin = direction;

        // We rotate the origin direction by the angular distribution angle
        G4ThreeVector orthoVector = direction.orthogonal();
        direction.rotate(angle, orthoVector);

        // We rotate a random angle along the radial direction
        Double_t randomAngle = G4UniformRand() * 2 * M_PI;
        direction.rotate(randomAngle, referenceOrigin);

        //       G4cout << "Angle  " << direction.angle( referenceOrigin ) << "
        //       should be = to " << angle << G4endl;

    } else if (angular_dist_type == g4_metadata_parameters::angular_dist_types::FLUX) {
        TVector3 v = p.GetMomentumDirection();

        v = v.Unit();

        direction.set(v.X(), v.Y(), v.Z());

    } else if (angular_dist_type == g4_metadata_parameters::angular_dist_types::BACK_TO_BACK) {
        // This should never crash. In TRestG4Metadata we have defined that if the
        // first source is backtoback we set it to isotropic
        TVector3 v = p.GetMomentumDirection();
        v = v.Unit();

        direction.set(-v.X(), -v.Y(), -v.Z());
    } else {
        G4cout << "WARNING: Generator angular distribution was not recognized. "
                  "Launching particle to (1,0,0)"
               << G4endl;
    }

    // setting particle direction
    fParticleGun->SetParticleMomentumDirection(direction);
}

void PrimaryGeneratorAction::SetParticleEnergy(Int_t n, TRestGeant4Particle p) {
    Double_t energy = 0;

    string energy_dist_type_name = (string)fRestGeant4Metadata->GetParticleSource(n)->GetEnergyDistType();
    energy_dist_type_name = g4_metadata_parameters::CleanString(energy_dist_type_name);

    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Debug) {
        G4cout << "DEBUG: Energy distribution: " << energy_dist_type_name << endl;
    }

    g4_metadata_parameters::energy_dist_types energy_dist_type;
    if (g4_metadata_parameters::energy_dist_types_map.count(energy_dist_type_name)) {
        energy_dist_type = g4_metadata_parameters::energy_dist_types_map[energy_dist_type_name];
    } else {
        // if we get here it means the parameter is not valid, we can either assign a default value or stop
        // execution default value in this case is 1 keV
        G4cout << "Invalid energy distribution (" + energy_dist_type_name + ") valid values are: ";
        for (const auto& pair : g4_metadata_parameters::energy_dist_types_map) {
            G4cout << pair.first << ", ";
        }
        G4cout << std::endl;
        G4cout << "WARNING! Energy distribution type was not recognized. Setting "
                  "energy to 1keV"
               << G4endl;
        energy = 1 * keV;
        // maybe it would be better to set type as mono and change energy of particle, or define the default
        // energy there throw "Invalid energy distribution type";
    }

    if (energy_dist_type == g4_metadata_parameters::energy_dist_types::MONO) {
        energy = p.GetEnergy() * keV;
    } else if (energy_dist_type == g4_metadata_parameters::energy_dist_types::FLAT) {
        TVector2 enRange = fRestGeant4Metadata->GetParticleSource(n)->GetEnergyRange();
        energy = ((enRange.Y() - enRange.X()) * G4UniformRand() + enRange.X()) * keV;
    } else if (energy_dist_type == g4_metadata_parameters::energy_dist_types::LOG) {
        TVector2 enRange = fRestGeant4Metadata->GetParticleSource(n)->GetEnergyRange();
        auto max_energy = enRange.Y() * keV;
        auto min_energy = enRange.X() * keV;
        energy = exp((log(max_energy) - log(min_energy)) * G4UniformRand() + log(min_energy));

    } else if (energy_dist_type == g4_metadata_parameters::energy_dist_types::TH1D) {
        Double_t value = G4UniformRand() * fSpectrumIntegral;
        Double_t sum = 0;
        Double_t deltaEnergy =
            fPrimaryEnergySpectrum.GetBinCenter(2) - fPrimaryEnergySpectrum.GetBinCenter(1);
        for (int bin = startEnergyBin; bin <= endEnergyBin; bin++) {
            sum += fPrimaryEnergySpectrum.GetBinContent(bin);

            if (sum > value) {
                energy = energyFactor *
                         (Double_t)(fPrimaryEnergySpectrum.GetBinCenter(bin) +
                                    deltaEnergy * (0.5 - G4UniformRand())) *
                         keV;
                break;
            }
        }
    } else {
        G4cout << "WARNING! Energy distribution type was not recognized. Setting "
                  "energy to 1keV"
               << G4endl;
        energy = 1 * keV;
    }

    string angular_dist_type_name = (string)fRestGeant4Metadata->GetParticleSource(n)->GetAngularDistType();
    angular_dist_type_name = g4_metadata_parameters::CleanString(angular_dist_type_name);
    g4_metadata_parameters::angular_dist_types angular_dist_type;
    if (g4_metadata_parameters::angular_dist_types_map.count(angular_dist_type_name)) {
        angular_dist_type = g4_metadata_parameters::angular_dist_types_map[angular_dist_type_name];
        if (n > 0 && angular_dist_type == g4_metadata_parameters::angular_dist_types::BACK_TO_BACK)
            energy = lastEnergy;
    }

    if (n == 0) lastEnergy = energy;
    fParticleGun->SetParticleEnergy(energy);

    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Debug)
        G4cout << "DEBUG: Particle energy: " << energy / keV << " keV" << endl;
}

void PrimaryGeneratorAction::SetParticlePosition() {
    double x = 0, y = 0, z = 0;
    string generator_type_name = (string)fRestGeant4Metadata->GetGeneratorType();
    g4_metadata_parameters::generator_types generator_type =
        g4_metadata_parameters::generator_types_map[g4_metadata_parameters::CleanString(generator_type_name)];

    string generator_shape_name = (string)fRestGeant4Metadata->GetGeneratorShape();
    g4_metadata_parameters::generator_shapes generator_shape =
        g4_metadata_parameters::generator_shapes_map[g4_metadata_parameters::CleanString(
            generator_shape_name)];

    while (true) {
        if (generator_type == g4_metadata_parameters::generator_types::POINT) {
            GenPositionOnPoint(x, y, z);
        } else if (generator_type == g4_metadata_parameters::generator_types::SURFACE) {
            if (generator_shape == g4_metadata_parameters::generator_shapes::GDML) {
                GenPositionOnGDMLSurface(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::BOX) {
                GenPositionOnBoxSurface(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::CYLINDER) {
                GenPositionOnCylinderSurface(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::SPHERE) {
                GenPositionOnSphereSurface(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::PLATE) {
                GenPositionOnPlate(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::WALL) {
                GenPositionOnWall(x, y, z);
            }
        } else if (generator_type == g4_metadata_parameters::generator_types::VOLUME) {
            if (generator_shape == g4_metadata_parameters::generator_shapes::GDML) {
                GenPositionOnGDMLVolume(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::BOX) {
                GenPositionOnBoxVolume(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::CYLINDER) {
                GenPositionOnCylinderVolume(x, y, z);
            } else if (generator_shape == g4_metadata_parameters::generator_shapes::SPHERE) {
                GenPositionOnSphereVolume(x, y, z);
            }
        } else {
            G4cout << "WARNING! Generator type \"" << fRestGeant4Metadata->GetGeneratorType()
                   << "\" was not recognized. Launching particle "
                      "from origin (0,0,0)"
                   << G4endl;
        }

        // use the density funciton. If the density is small, then val2 is small, we are more
        // likely to regenerate the particle position
        if (fGeneratorSpatialDensityFunction) {
            double val1 = G4UniformRand();
            double val2 = fGeneratorSpatialDensityFunction->Eval(x, y, z);
            if (val2 > 1) {
                G4cout << "error! Generator density function > 1 at position (" << x << ", " << y << ", " << z
                     << "), check your definition!" << endl;
                abort();
            }
            if (val1 > val2) {
                continue;
            }
        }
        break;
    }

    // storing the direction in TRestG4Event class
    TVector3 eventPosition(x, y, z);
    // restG4Event->SetPrimaryEventOrigin(eventPosition);

    // setting particle position
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
}

////_____________________________________________________________________________
// void PrimaryGeneratorAction::SetParticlePosition(int n) {
//    // Storing particle's position to that retrieved from TRestG4Particle
//    TVector3 pos = p.GetOrigin();
//    restG4Event->SetPrimaryEventOrigin(pos);
//
//    if (fRestGeant4Metadata->GetVerboseLevel() >= REST_Debug) {
//        G4cout << "DEBUG: Event origin: "
//             << "(" << restG4Event->GetPrimaryEventOrigin().X() << ", "
//             << restG4Event->GetPrimaryEventOrigin().Y() << ", " << restG4Event->GetPrimaryEventOrigin().Z()
//             << ")"
//             << " mm" << endl;
//    }
//
//    // Setting particle position
//    fParticleGun->SetParticlePosition(G4ThreeVector(pos.X(), pos.Y(), pos.Z()));
//}

G4ThreeVector PrimaryGeneratorAction::GetIsotropicVector() {
    G4double a, b, c;
    G4double n;
    do {
        a = (G4UniformRand() - 0.5) / 0.5;
        b = (G4UniformRand() - 0.5) / 0.5;
        c = (G4UniformRand() - 0.5) / 0.5;
        n = a * a + b * b + c * c;
    } while (n > 1 || n == 0.0);

    n = std::sqrt(n);
    a /= n;
    b /= n;
    c /= n;
    return G4ThreeVector(a, b, c);
}

//
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
    double xMin = fDetector->GetBoundingX_min();
    double xMax = fDetector->GetBoundingX_max();
    double yMin = fDetector->GetBoundingY_min();
    double yMax = fDetector->GetBoundingY_max();
    double zMin = fDetector->GetBoundingZ_min();
    double zMax = fDetector->GetBoundingZ_max();

    do {
        x = xMin + (xMax - xMin) * G4UniformRand();
        y = yMin + (yMax - yMin) * G4UniformRand();
        z = zMin + (zMax - zMin) * G4UniformRand();
    } while (fDetector->GetGeneratorSolid()->Inside(G4ThreeVector(x, y, z)) != kInside);

    x = x + fDetector->GetGeneratorTranslation().x();
    y = y + fDetector->GetGeneratorTranslation().y();
    z = z + fDetector->GetGeneratorTranslation().z();
}
void PrimaryGeneratorAction::GenPositionOnGDMLSurface(double& x, double& y, double& z) {
    // TODO there is a problem, probably with G4 function GetPointOnSurface
    // It produces a point on the surface but it is not uniformly distributed
    // May be it is just an OPENGL drawing issue?

    G4ThreeVector position = fDetector->GetGeneratorSolid()->GetPointOnSurface();

    x = position.x();
    y = position.y();
    z = position.z();

    x = x + fDetector->GetGeneratorTranslation().x();
    y = y + fDetector->GetGeneratorTranslation().y();
    z = z + fDetector->GetGeneratorTranslation().z();
}
void PrimaryGeneratorAction::GenPositionOnBoxVolume(double& x, double& y, double& z) {
    Double_t sidex = fRestGeant4Metadata->GetGeneratorSize().X();
    Double_t sidey = fRestGeant4Metadata->GetGeneratorSize().Y();
    Double_t sidez = fRestGeant4Metadata->GetGeneratorSize().Z();

    x = sidex * (G4UniformRand() - 0.5);
    y = sidey * (G4UniformRand() - 0.5);
    z = sidez * (G4UniformRand() - 0.5);

    G4ThreeVector rndPos = G4ThreeVector(x, y, z);

    TVector3 rotaxis = fRestGeant4Metadata->GetGeneratorRotationAxis();
    G4ThreeVector rotaxisG4 = G4ThreeVector(rotaxis.x(), rotaxis.y(), rotaxis.z());
    rndPos.rotate(rotaxisG4, fRestGeant4Metadata->GetGeneratorRotationDegree());

    TVector3 center = fRestGeant4Metadata->GetGeneratorPosition();

    x = rndPos.x() + center.X();
    y = rndPos.y() + center.Y();
    z = rndPos.z() + center.Z();
}
void PrimaryGeneratorAction::GenPositionOnBoxSurface(double& x, double& y, double& z) {
    G4cout << __PRETTY_FUNCTION__ << ": not implemented!" << endl;
    abort();
}
void PrimaryGeneratorAction::GenPositionOnSphereVolume(double& x, double& y, double& z) {
    G4cout << __PRETTY_FUNCTION__ << ": not implemented!" << endl;
    abort();
}
void PrimaryGeneratorAction::GenPositionOnSphereSurface(double& x, double& y, double& z) {
    G4ThreeVector rndPos = GetIsotropicVector();

    Double_t radius = fRestGeant4Metadata->GetGeneratorSize().X();

    TVector3 center = fRestGeant4Metadata->GetGeneratorPosition();

    x = radius * rndPos.x() + center.X();
    y = radius * rndPos.y() + center.Y();
    z = radius * rndPos.z() + center.Z();
}
void PrimaryGeneratorAction::GenPositionOnCylinderVolume(double& x, double& y, double& z) {
    G4cout << __PRETTY_FUNCTION__ << ": not implemented!" << endl;
    abort();
}
void PrimaryGeneratorAction::GenPositionOnCylinderSurface(double& x, double& y, double& z) {
    Double_t angle = 2 * M_PI * G4UniformRand();

    Double_t radius = fRestGeant4Metadata->GetGeneratorSize().x();
    Double_t length = fRestGeant4Metadata->GetGeneratorSize().y();

    x = radius * cos(angle);
    y = radius * sin(angle);
    z = length * (G4UniformRand() - 0.5);

    G4ThreeVector rndPos = G4ThreeVector(x, y, z);

    TVector3 rotaxis = fRestGeant4Metadata->GetGeneratorRotationAxis();
    G4ThreeVector rotaxisG4 = G4ThreeVector(rotaxis.x(), rotaxis.y(), rotaxis.z());
    rndPos.rotate(rotaxisG4, fRestGeant4Metadata->GetGeneratorRotationDegree());

    TVector3 center = fRestGeant4Metadata->GetGeneratorPosition();

    x = rndPos.x() + center.X();
    y = rndPos.y() + center.Y();
    z = rndPos.z() + center.Z();
}
void PrimaryGeneratorAction::GenPositionOnPoint(double& x, double& y, double& z) {
    TVector3 position = fRestGeant4Metadata->GetGeneratorPosition();

    x = position.X();
    y = position.Y();
    z = position.Z();
}
void PrimaryGeneratorAction::GenPositionOnWall(double& x, double& y, double& z) {
    Double_t sidex = fRestGeant4Metadata->GetGeneratorSize().X();
    Double_t sidey = fRestGeant4Metadata->GetGeneratorSize().Y();

    x = sidex * (G4UniformRand() - 0.5);
    y = sidey * (G4UniformRand() - 0.5);

    G4ThreeVector rndPos = G4ThreeVector(x, y, 0);

    TVector3 rotaxis = fRestGeant4Metadata->GetGeneratorRotationAxis();
    G4ThreeVector rotaxisG4 = G4ThreeVector(rotaxis.x(), rotaxis.y(), rotaxis.z());
    rndPos.rotate(rotaxisG4, fRestGeant4Metadata->GetGeneratorRotationDegree());

    TVector3 center = fRestGeant4Metadata->GetGeneratorPosition();

    x = rndPos.x() + center.X();
    y = rndPos.y() + center.Y();
    z = rndPos.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnPlate(double& x, double& y, double& z) {
    Double_t radius = fRestGeant4Metadata->GetGeneratorSize().X();

    do {
        x = 2 * radius * (G4UniformRand() - 0.5);
        y = 2 * radius * (G4UniformRand() - 0.5);
        //       G4cout << "x : " << x << " y : " << y << endl;
    } while (x * x + y * y > radius * radius);

    G4ThreeVector rndPos = G4ThreeVector(x, y, 0);

    TVector3 rotaxis = fRestGeant4Metadata->GetGeneratorRotationAxis();
    G4ThreeVector rotaxisG4 = G4ThreeVector(rotaxis.x(), rotaxis.y(), rotaxis.z());
    rndPos.rotate(rotaxisG4, fRestGeant4Metadata->GetGeneratorRotationDegree());

    TVector3 center = fRestGeant4Metadata->GetGeneratorPosition();

    x = rndPos.x() + center.X();
    y = rndPos.y() + center.Y();
    z = rndPos.z() + center.Z();
}
