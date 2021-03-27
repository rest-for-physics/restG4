#include "DetectorConstruction.hh"

#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4UserLimits.hh>

#include "G4FieldManager.hh"
#include "G4IonTable.hh"
#include "G4Isotope.hh"
#include "G4MagneticField.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"

extern TRestGeant4Event* restG4Event;
extern TRestGeant4Metadata* restG4Metadata;

//_____________________________________________________________________________
DetectorConstruction::DetectorConstruction() {
    G4cout << "Detector Construction" << G4endl;

    parser = new G4GDMLParser();
}
//_____________________________________________________________________________
DetectorConstruction::~DetectorConstruction() { delete parser; }
//_____________________________________________________________________________
G4VPhysicalVolume* DetectorConstruction::Construct() {
    cout << "Isotope table " << endl;
    cout << *(G4Isotope::GetIsotopeTable()) << endl;

    G4cout << "Producing geometry" << G4endl;

    // Reading the geometry
    TString geometryFile = restG4Metadata->Get_GDML_Filename();

    // char originDirectory[255];
    // sprintf( originDirectory, "%s", get_current_dir_name() );

    // char buffer[255];
    // sprintf( buffer, "%s", (char *) restG4Metadata->GetGeometryPath().Data() );
    // chdir( buffer );
    char originDirectory[256];
    sprintf(originDirectory, "%s", getenv("PWD"));
    auto pathandname = TRestTools::SeparatePathAndName((string)restG4Metadata->Get_GDML_Filename());
    chdir(pathandname.first.c_str());

    parser->Read(pathandname.second, false);

    chdir(originDirectory);

    G4VPhysicalVolume* W = parser->GetWorldVolume();

    // TODO : Take the name of the sensitive volume and use it here to define its
    // StepSize
    string sensibleVolumeName = (string)restG4Metadata->GetSensitiveVolume();
    bool foundSensitiveVolume = false;
    G4PhysicalVolumeStore* physicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
    cout << "Physical volumes (Geant4)" << endl;
    for (const auto& physicalVolume : *physicalVolumeStore) {
        cout << "\t" << (string)physicalVolume->GetName() << endl;
        if ((string)physicalVolume->GetName() == sensibleVolumeName) {
            foundSensitiveVolume = true;
            break;
        }
    }

    if (!foundSensitiveVolume) {
        // If using GDML assembly, the names of G4VPhysicalVolume and GDML physical volume do not match, and I
        // have found no way to recover the GDML physical volume name using this parser (if you find one,
        // PLEASE update / write an issue / pull request)
        //
        // one workaround is using the logical volume name as an identifier of the physical volume, this works
        // if the logical volume is only used in one physical volume

        // do an attempt to match logical volume to physical volume before exiting

        map<string, string> physicalToLogicalMap;
        for (const auto& physicalVolume : *physicalVolumeStore) {
            physicalToLogicalMap[physicalVolume->GetName()] = physicalVolume->GetLogicalVolume()->GetName();
        }

        int count = 0;
        string sensitivePhysicalVolumeName;
        for (auto const& kv : physicalToLogicalMap) {
            if (kv.second == sensibleVolumeName) {
                count++;
                sensitivePhysicalVolumeName = kv.first;
            }
        }

        if (count == 1) {
            // this means the logical volume only appeared once, so we can safely get the physical volume name
            cout << "WARNING: Getting the physical volume name ('" << sensitivePhysicalVolumeName
                 << "') from LOGICAL VOLUME name: " << sensibleVolumeName << endl;
            cout << "WARNING: Changing RestG4Metadata Sensitive Volume from '" << sensibleVolumeName
                 << "' to: '" << sensitivePhysicalVolumeName << "'" << endl;
            restG4Metadata->SetSensitiveVolume((TString)sensitivePhysicalVolumeName);
            sensibleVolumeName = (string)restG4Metadata->GetSensitiveVolume();
        } else {
            G4cout << "RESTG4 error. Sensitive volume '" << sensibleVolumeName
                   << "' does not exist in geometry!!" << G4endl;
            G4cout << "RESTG4 error. Please, review geometry! Press a key to crash!!" << G4endl;
            getchar();
            // We need to produce a clean exit at this point
            exit(1);
        }
    }
    G4VPhysicalVolume* _vol = GetPhysicalVolume(sensibleVolumeName);

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
        cout << "ERROR : Logical volume for sensitive \"" << sensibleVolumeName << "\" not found!" << endl;
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
        TString activeVolumeName =
            restG4Metadata->GetActiveVolumeName(id);  // this is the name of the physical volume

        G4LogicalVolumeStore* logicalVolumeStore = G4LogicalVolumeStore::GetInstance();
        for (const auto& logicalVolume : *logicalVolumeStore) {
            cout << "\t" << logicalVolume->GetName() << endl;
            if (restG4Metadata->GetMaxStepSize(activeVolumeName) > 0) {
                G4cout << "Setting maxStepSize = " << restG4Metadata->GetMaxStepSize(activeVolumeName)
                       << "mm for volume : " << activeVolumeName << G4endl;
                logicalVolume->SetUserLimits(
                    new G4UserLimits(restG4Metadata->GetMaxStepSize(activeVolumeName) * mm));
            }

            cout << "Activating volume : " << activeVolumeName << endl;
            restG4Event->AddActiveVolume((string)activeVolumeName);
            if (logicalVolume == nullptr) {
                cout << "Error: DetectorConstruction. Logical volume " << activeVolumeName
                     << " is not defined in the geometry" << endl;
                exit(1);
            }
        }
    }

    cout << "Detector constructed : " << W << endl;

    return W;
}

G4VPhysicalVolume* DetectorConstruction::GetPhysicalVolume(G4String physicalVolumeName) {
    G4PhysicalVolumeStore* physicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
    for (const auto& physicalVolume : *physicalVolumeStore) {
        cout << "\t" << physicalVolume->GetName() << endl;
        if (physicalVolume->GetName() == physicalVolumeName) return physicalVolume;
    }

    return nullptr;
}
