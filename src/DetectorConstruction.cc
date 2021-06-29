
#include "DetectorConstruction.hh"

#include <G4FieldManager.hh>
#include <G4IonTable.hh>
#include <G4Isotope.hh>
#include <G4MagneticField.hh>
#include <G4Material.hh>
#include <G4SystemOfUnits.hh>
#include <G4UniformMagField.hh>
#include <G4UserLimits.hh>

extern TRestGeant4Event* restG4Event;
extern TRestGeant4Metadata* restG4Metadata;

DetectorConstruction::DetectorConstruction() {
    G4cout << "Detector Construction" << G4endl;
    parser = new G4GDMLParser();
}

DetectorConstruction::~DetectorConstruction() { delete parser; }

G4VPhysicalVolume* DetectorConstruction::Construct() {
    cout << "Isotope table " << endl;
    cout << *(G4Isotope::GetIsotopeTable()) << endl;

    G4cout << "Producing geometry" << G4endl;

    // Reading the geometry
    TString geometryFile = restG4Metadata->Get_GDML_Filename();

    char originDirectory[256];
    sprintf(originDirectory, "%s", getenv("PWD"));
    auto separatePathAndName = TRestTools::SeparatePathAndName((string)restG4Metadata->Get_GDML_Filename());
    chdir(separatePathAndName.first.c_str());

    parser->Read(separatePathAndName.second, false);

    chdir(originDirectory);

    G4VPhysicalVolume* W = parser->GetWorldVolume();

    // TODO : Take the name of the sensitive volume and use it here to define its
    // StepSize
    string SensVol = (string)restG4Metadata->GetSensitiveVolume();
    G4VPhysicalVolume* _vol = GetPhysicalVolume(SensVol);
    if (!_vol) {
        G4cout << "RESTG4 error. Sensitive volume  " << SensVol << " does not exist in geomtry!!" << G4endl;
        G4cout << "RESTG4 error. Please, review geometry! Presh a key to crash!!" << G4endl;
        getchar();
        // We need to produce a clean exit at this point
        exit(1);
    }

    Double_t mx = restG4Metadata->GetMagneticField().X() * tesla;
    Double_t my = restG4Metadata->GetMagneticField().Y() * tesla;
    Double_t mz = restG4Metadata->GetMagneticField().Z() * tesla;

    G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(mx, my, mz));
    G4FieldManager* localFieldMgr = new G4FieldManager(magField);
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);

    if (_vol != nullptr) {
        G4LogicalVolume* vol = _vol->GetLogicalVolume();
        // This method seems not available in my Geant4 version 10.4.2
        // In future Geant4 versions it seems possible to define field at particular volumes
        // vol->setFieldManager(localFieldMgr, true);
        G4Material* mat = vol->GetMaterial();
        G4cout << "Sensitivity volume properties" << G4endl;
        G4cout << "==============" << G4endl;
        G4cout << "Sensitivity volume name : " << mat->GetName() << G4endl;
        G4cout << "Sensitivity volume temperature : " << mat->GetTemperature() << G4endl;
        G4cout << "Sensitivity volume density : " << mat->GetDensity() / (g / cm3) << " g/cm3" << G4endl;
    } else {
        cout << "ERROR : Logical volume for sensitive \"" << SensVol << "\" not found!" << endl;
    }

    // Getting generation volume
    string GenVol = (string)restG4Metadata->GetGeneratedFrom();
    cout << "Generated from volume : " << GenVol << endl;
    string type = (string)restG4Metadata->GetGeneratorType();
    cout << "Generator type : " << type << endl;

    // TODO if we do not find the volume given in the config inside the geometry
    // we should RETURN error
    if (type == "volume" && GenVol != "Not defined") {
        G4VPhysicalVolume* pVol = GetPhysicalVolume(GenVol);
        if (pVol == nullptr) {
            cout << "ERROR : The generator volume was not found in the geometry" << endl;
            exit(1);
            return W;
        }

        fGeneratorTranslation = pVol->GetTranslation();

        // We set in TRestGeant4Metadata the center of the generator. If it is a point
        // we just want the value from the config file.
        // TODO : make this kind of keyword comparisons case insensitive?
        if (type == "surface" || type == "volume") {
            restG4Metadata->SetGeneratorPosition(fGeneratorTranslation.x(), fGeneratorTranslation.y(),
                                                 fGeneratorTranslation.z());
        }

        generatorSolid = pVol->GetLogicalVolume()->GetSolid();

        // while ( fDetector->GetGeneratorSolid()->Inside( G4ThreeVector( x, y, z) )
        // != kInside );

        // We estimate the maximum distance of our volume
        // The code below returns a value bigger than expected
        // If we try with a cylinder the maximum distance should be sqrt(R*R+L*L)
        // But the value returned by this is bigger TODO check this
        boundBox_xMax = -1.e30;
        boundBox_yMax = -1.e30;
        boundBox_zMax = -1.e30;
        boundBox_xMin = 1.e30;
        boundBox_yMin = 1.e30;
        boundBox_zMin = 1.e30;
        if (type == "volume") {
            cout << "Optimizing REST volume generation (Please wait. This might take "
                    "few minutes depending on geometry complexity) "
                 << flush;

            for (int n = 0; n < 100000; n++) {
                G4ThreeVector point = generatorSolid->GetPointOnSurface();

                if (point.x() > boundBox_xMax) boundBox_xMax = point.x();
                if (point.y() > boundBox_yMax) boundBox_yMax = point.y();
                if (point.z() > boundBox_zMax) boundBox_zMax = point.z();

                if (point.x() < boundBox_xMin) boundBox_xMin = point.x();
                if (point.y() < boundBox_yMin) boundBox_yMin = point.y();
                if (point.z() < boundBox_zMin) boundBox_zMin = point.z();
            }

            boundBox_xMin = boundBox_xMin * 1.1;
            boundBox_xMax = boundBox_xMax * 1.1;

            boundBox_yMin = boundBox_yMin * 1.1;
            boundBox_yMax = boundBox_yMax * 1.1;

            boundBox_zMin = boundBox_zMin * 1.1;
            boundBox_zMax = boundBox_zMax * 1.1;
        }
    }

    for (int id = 0; id < restG4Metadata->GetNumberOfActiveVolumes(); id++) {
        TString actVolName = restG4Metadata->GetActiveVolumeName(id);
        G4VPhysicalVolume* pVol = GetPhysicalVolume((G4String)actVolName);
        if (pVol != nullptr) {
            G4LogicalVolume* lVol = pVol->GetLogicalVolume();
            if (restG4Metadata->GetMaxStepSize(actVolName) > 0) {
                G4cout << "Setting maxStepSize = " << restG4Metadata->GetMaxStepSize(actVolName)
                       << "mm for volume : " << actVolName << G4endl;
                lVol->SetUserLimits(new G4UserLimits(restG4Metadata->GetMaxStepSize(actVolName) * mm));
            }
        }

        cout << "Activating volume : " << actVolName << endl;
        restG4Event->AddActiveVolume((string)actVolName);
        if (pVol == nullptr) {
            cout << "DetectorConstruction. Volume " << actVolName << " is not defined in the geometry"
                 << endl;
            exit(1);
        }
    }

    cout << "Detector constructed : " << W << endl;

    return W;
}

G4VPhysicalVolume* DetectorConstruction::GetPhysicalVolume(G4String physVolName) {
    G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
    std::vector<G4VPhysicalVolume*>::const_iterator physVol;
    for (physVol = physVolStore->begin(); physVol != physVolStore->end(); physVol++) {
        if ((*physVol)->GetName() == physVolName) {
            return *physVol;
        }
    }

    return nullptr;
}
