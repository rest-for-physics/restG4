
#include "DetectorConstruction.h"

#include <SensitiveDetector.h>
#include <TRestGeant4GeometryInfo.h>

#include <G4FieldManager.hh>
#include <G4IonTable.hh>
#include <G4Isotope.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4MagneticField.hh>
#include <G4Material.hh>
#include <G4RunManager.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UniformMagField.hh>
#include <G4UserLimits.hh>
//optical surfaces
#include <G4LogicalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>
#include <filesystem>

#include "SimulationManager.h"

using namespace std;
using namespace TRestGeant4PrimaryGeneratorTypes;

DetectorConstruction::DetectorConstruction(SimulationManager* simulationManager)
    : fSimulationManager(simulationManager) {
    G4cout << "Detector Construction" << G4endl;
    fGdmlParser = new G4GDMLParser();
}

DetectorConstruction::~DetectorConstruction() { delete fGdmlParser; }

G4VPhysicalVolume* DetectorConstruction::Construct() {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();

    cout << "Isotope table " << endl;
    cout << *(G4Isotope::GetIsotopeTable()) << endl;

    G4cout << "Producing geometry" << G4endl;

    // Reading the geometry
    TString geometryFile = restG4Metadata->GetGdmlFilename();

    const auto startingPath = filesystem::current_path();

    const auto [gdmlPath, gdmlToRead] =
        TRestTools::SeparatePathAndName((string)restG4Metadata->GetGdmlFilename());
    filesystem::current_path(gdmlPath);

    G4cout << "gdmlToRead: " << gdmlToRead << G4endl;

    fGdmlParser->Read(gdmlToRead, false);
    G4VPhysicalVolume* worldVolume = fGdmlParser->GetWorldVolume();

    const auto worldSolid = dynamic_cast<G4Box*>(worldVolume->GetLogicalVolume()->GetSolid());
    restG4Metadata->fGeant4PrimaryGeneratorInfo.fSpatialGeneratorWorldSize = {
        worldSolid->GetXHalfLength(), worldSolid->GetYHalfLength(), worldSolid->GetZHalfLength()};

    restG4Metadata->fGeant4GeometryInfo.InitializeOnDetectorConstruction(gdmlToRead, worldVolume);
    restG4Metadata->ReadDetector();
    restG4Metadata->PrintMetadata();  // now we have detector info

    const auto& geometryInfo = restG4Metadata->GetGeant4GeometryInfo();
    geometryInfo.Print();

    // do some checks
    {
        // Check all physical volume names are unique
        G4PhysicalVolumeStore* physicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
        set<string> physicalVolumeNames;
        vector<G4VPhysicalVolume*>::const_iterator physicalVolume;
        for (physicalVolume = physicalVolumeStore->begin(); physicalVolume != physicalVolumeStore->end();
             physicalVolume++) {
            auto name = (*physicalVolume)->GetName();
            if (physicalVolumeNames.count(name)) {
                cerr << "ERROR: physical volume name " << name
                     << " is not unique. Please double check your geometry files. Be mindful of especial "
                        "characters such as '0x'"
                     << endl;
                exit(1);
            }
            physicalVolumeNames.insert(name);
        }

        // Check all logical volume names are unique
        G4LogicalVolumeStore* logicalVolumeStore = G4LogicalVolumeStore::GetInstance();
        set<string> logicalVolumeNames;
        vector<G4LogicalVolume*>::const_iterator logicalVolume;
        for (logicalVolume = logicalVolumeStore->begin(); logicalVolume != logicalVolumeStore->end();
             logicalVolume++) {
            auto name = (*logicalVolume)->GetName();
            if (logicalVolumeNames.count(name)) {
                cerr << "ERROR: logical volume name " << name
                     << " is not unique. Please double check your geometry files. Be mindful of especial "
                        "characters such as '0x'"
                     << endl;
                exit(1);
            }
            logicalVolumeNames.insert(name);
        }
    }
    filesystem::current_path(startingPath);

    auto sensitiveVolume = (string)restG4Metadata->GetSensitiveVolume();
    G4VPhysicalVolume* sensitivePhysicalVolume = GetPhysicalVolume(sensitiveVolume);
    if (sensitivePhysicalVolume == nullptr) {
        // sensitive volume was not found, perhaps the user specified a logical volume
        auto physicalVolumes = geometryInfo.GetAllPhysicalVolumesFromLogical(sensitiveVolume);
        if (physicalVolumes.size() == 1) {
            restG4Metadata->InsertSensitiveVolume(
                geometryInfo.GetAlternativeNameFromGeant4PhysicalName(physicalVolumes[0]));
            sensitiveVolume = (string)restG4Metadata->GetSensitiveVolume();
            sensitivePhysicalVolume = GetPhysicalVolume(sensitiveVolume);
        }
    }

    if (sensitivePhysicalVolume == nullptr) {
        cerr << "ERROR: Sensitive volume '" << sensitiveVolume << "' not found" << endl;
        exit(1);
    }

    Double_t mx = restG4Metadata->GetMagneticField().X() * tesla;
    Double_t my = restG4Metadata->GetMagneticField().Y() * tesla;
    Double_t mz = restG4Metadata->GetMagneticField().Z() * tesla;

    G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(mx, my, mz));
    // G4FieldManager* localFieldMgr = new G4FieldManager(magField);
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);

    G4LogicalVolume* volume = sensitivePhysicalVolume->GetLogicalVolume();
    G4Material* material = volume->GetMaterial();
    G4cout << "Sensitive volume properties:" << G4endl;
    G4cout << "\t- Material: " << material->GetName() << G4endl;
    G4cout << "\t- Temperature: " << material->GetTemperature() << " K" << G4endl;
    G4cout << "\t- Density: " << material->GetDensity() / (g / cm3) << " g/cm3" << G4endl;

    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();
    // Getting generation volume
    const auto fromVolume = primaryGeneratorInfo.GetSpatialGeneratorFrom();
    if (fromVolume != "NO_SUCH_PARA") {
        cout << "Generated from volume: " << primaryGeneratorInfo.GetSpatialGeneratorFrom() << endl;
    }
    cout << "Generator type: " << primaryGeneratorInfo.GetSpatialGeneratorType() << endl;

    const auto spatialGeneratorTypeEnum =
        StringToSpatialGeneratorTypes(primaryGeneratorInfo.GetSpatialGeneratorType().Data());

    if (spatialGeneratorTypeEnum == TRestGeant4PrimaryGeneratorTypes::SpatialGeneratorTypes::VOLUME &&
        primaryGeneratorInfo.GetSpatialGeneratorFrom() != "Not defined") {
        G4VPhysicalVolume* gdmlPhysicalVolume =
            GetPhysicalVolume(primaryGeneratorInfo.GetSpatialGeneratorFrom().Data());
        if (gdmlPhysicalVolume == nullptr) {
            // perhaps the user selected a logical volume instead
            auto physicalVolumes = geometryInfo.GetAllPhysicalVolumesFromLogical(
                primaryGeneratorInfo.GetSpatialGeneratorFrom().Data());
            if (physicalVolumes.size() == 1) {
                gdmlPhysicalVolume = GetPhysicalVolume(physicalVolumes[0].Data());
                cout << "Generator volume '" << primaryGeneratorInfo.GetSpatialGeneratorFrom()
                     << "' was not found in the geometry. Using the physical volume '" << physicalVolumes[0]
                     << "' instead, which was obtained from logical volume '"
                     << primaryGeneratorInfo.GetSpatialGeneratorFrom() << "'" << endl;
            }
        }
        if (gdmlPhysicalVolume == nullptr) {
            cerr << "ERROR: The generator volume '" << primaryGeneratorInfo.GetSpatialGeneratorFrom()
                 << "' was not found in the geometry" << endl;
            exit(1);
        }

        fGeneratorTranslation = gdmlPhysicalVolume->GetTranslation();
        if (spatialGeneratorTypeEnum == TRestGeant4PrimaryGeneratorTypes::SpatialGeneratorTypes::SURFACE ||
            spatialGeneratorTypeEnum == TRestGeant4PrimaryGeneratorTypes::SpatialGeneratorTypes::VOLUME) {
            restG4Metadata->fGeant4PrimaryGeneratorInfo.fSpatialGeneratorPosition = {
                fGeneratorTranslation.x(), fGeneratorTranslation.y(), fGeneratorTranslation.z()};
        }

        fGeneratorSolid = gdmlPhysicalVolume->GetLogicalVolume()->GetSolid();

        fBoundBoxXMax = -1.e30;
        fBoundBoxYMax = -1.e30;
        fBoundBoxZMax = -1.e30;
        fBoundBoxXMin = 1.e30;
        fBoundBoxYMin = 1.e30;
        fBoundBoxZMin = 1.e30;
        if (spatialGeneratorTypeEnum == TRestGeant4PrimaryGeneratorTypes::SpatialGeneratorTypes::VOLUME) {
            cout << "Optimizing REST volume generation (Please wait. This might take "
                    "few minutes depending on geometry complexity) "
                 << flush;

            for (int n = 0; n < 100000; n++) {
                G4ThreeVector point = fGeneratorSolid->GetPointOnSurface();

                if (point.x() > fBoundBoxXMax) fBoundBoxXMax = point.x();
                if (point.y() > fBoundBoxYMax) fBoundBoxYMax = point.y();
                if (point.z() > fBoundBoxZMax) fBoundBoxZMax = point.z();

                if (point.x() < fBoundBoxXMin) fBoundBoxXMin = point.x();
                if (point.y() < fBoundBoxYMin) fBoundBoxYMin = point.y();
                if (point.z() < fBoundBoxZMin) fBoundBoxZMin = point.z();
            }

            fBoundBoxXMin = fBoundBoxXMin * 1.1;
            fBoundBoxXMax = fBoundBoxXMax * 1.1;

            fBoundBoxYMin = fBoundBoxYMin * 1.1;
            fBoundBoxYMax = fBoundBoxYMax * 1.1;

            fBoundBoxZMin = fBoundBoxZMin * 1.1;
            fBoundBoxZMax = fBoundBoxZMax * 1.1;
        }
    }

    for (unsigned int id = 0; id < restG4Metadata->GetNumberOfActiveVolumes(); id++) {
        TString activeVolumeName = restG4Metadata->GetActiveVolumeName(id);
        G4VPhysicalVolume* physicalVolume = GetPhysicalVolume((G4String)activeVolumeName);
        if (physicalVolume != nullptr) {
            G4LogicalVolume* logicalVolume = physicalVolume->GetLogicalVolume();
            if (restG4Metadata->GetMaxStepSize(activeVolumeName) > 0) {
                RESTInfo << "Setting maxStepSize of " << restG4Metadata->GetMaxStepSize(activeVolumeName) * mm
                         << "mm for volume '" << activeVolumeName << "'" << RESTendl;
                logicalVolume->SetUserLimits(
                    new G4UserLimits(restG4Metadata->GetMaxStepSize(activeVolumeName) * mm));
            }
        } else {
            cerr << "DetectorConstruction::Construct - Volume '" << activeVolumeName
                 << "' is not defined in the geometry" << endl;
            exit(1);
        }
    }

    return worldVolume;
}

G4VPhysicalVolume* DetectorConstruction::GetPhysicalVolume(const G4String& physicalVolumeName) const {
    G4PhysicalVolumeStore* physicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& geometryInfo = restG4Metadata->GetGeant4GeometryInfo();
    vector<G4VPhysicalVolume*>::const_iterator physicalVolume;
    for (physicalVolume = physicalVolumeStore->begin(); physicalVolume != physicalVolumeStore->end();
         physicalVolume++) {
        auto name = (*physicalVolume)->GetName();
        auto alternativeName = (G4String)geometryInfo.GetAlternativeNameFromGeant4PhysicalName(name);
        if (name == physicalVolumeName || alternativeName == physicalVolumeName) {
            return *physicalVolume;
        }
    }

    return nullptr;
}

void DetectorConstruction::ConstructSDandField() {
    const TRestGeant4Metadata& metadata = *fSimulationManager->GetRestMetadata();

    set<G4LogicalVolume*> logicalVolumesSelected;
    for (const auto& userSensitiveVolume : metadata.GetSensitiveVolumes()) {
        // Each sensitive detector declaration in the RML should correspond to at least one logical volume
        G4LogicalVolume* logicalVolume;
        // Check if user selected a Geant4 physical volume by name
        G4VPhysicalVolume* physicalVolume =
            G4PhysicalVolumeStore::GetInstance()->GetVolume(userSensitiveVolume.Data(), false);
        if (physicalVolume == nullptr) {
            const G4String geant4VolumeName =
                metadata.GetGeant4GeometryInfo()
                    .GetGeant4PhysicalNameFromAlternativeName(userSensitiveVolume.Data())
                    .Data();
            // Check if user selected a Geant4 physical volume by REST name (only relevant for assemblies)
            physicalVolume = G4PhysicalVolumeStore::GetInstance()->GetVolume(geant4VolumeName, false);
        }
        if (physicalVolume == nullptr) {
            // Check if user selected a Geant4 logical volume by name
            logicalVolume = G4LogicalVolumeStore::GetInstance()->GetVolume(userSensitiveVolume.Data(), false);
        } else {
            // Sensitive detector already found
            logicalVolume = physicalVolume->GetLogicalVolume();
        }
        if (logicalVolume != nullptr) {
            logicalVolumesSelected.insert(logicalVolume);
        } else {
            // Check if the user string matches a logical volume by expanding the input as a regex
            // This can have multiple hits
            const auto logicalVolumesMatchingExpression =
                metadata.GetGeant4GeometryInfo().GetAllLogicalVolumesMatchingExpression(userSensitiveVolume);
            if (logicalVolumesMatchingExpression.empty()) {
                cerr << "Detector construction error: could not find matching logical volume(s) for '"
                     << userSensitiveVolume << "'" << endl;
                exit(1);
            } else {
                for (const auto& logicalVolumeName : logicalVolumesMatchingExpression) {
                    logicalVolumesSelected.insert(
                        G4LogicalVolumeStore::GetInstance()->GetVolume(logicalVolumeName.Data(), false));
                }
            }
        }
    }

    G4SDManager* SDManager = G4SDManager::GetSDMpointer();

    for (G4LogicalVolume* logicalVolume : logicalVolumesSelected) {
        auto name = logicalVolume->GetName();
        G4VSensitiveDetector* sensitiveDetector = new SensitiveDetector(fSimulationManager, name);
        SDManager->AddNewDetector(sensitiveDetector);
        logicalVolume->SetSensitiveDetector(sensitiveDetector);

        auto region = new G4Region(name);
        logicalVolume->SetRegion(region);
    }
}

void TRestGeant4GeometryInfo::PopulateFromGeant4World(const G4VPhysicalVolume* world) {
    auto detector = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    TRestGeant4Metadata* restG4Metadata = detector->fSimulationManager->GetRestMetadata();

    const size_t n = int(world->GetLogicalVolume()->GetNoDaughters());
    for (size_t i = 0; i < n + 1; i++) {  // world is the + 1
        G4VPhysicalVolume* volume;
        if (i == n) {
            volume = const_cast<G4VPhysicalVolume*>(world);
        } else {
            volume = world->GetLogicalVolume()->GetDaughter(i);
        }
        TString namePhysical = (TString)volume->GetName();
        if (fGdmlNewPhysicalNames.size() > i) {
            // it has been filled
            fGeant4PhysicalNameToNewPhysicalNameMap[namePhysical] = fGdmlNewPhysicalNames[i];
        }
        TString physicalNewName = GetAlternativeNameFromGeant4PhysicalName(namePhysical);
        TString nameLogical = (TString)volume->GetLogicalVolume()->GetName();
        TString nameMaterial = (TString)volume->GetLogicalVolume()->GetMaterial()->GetName();
        auto position = volume->GetTranslation();

        fPhysicalToLogicalVolumeMap[physicalNewName] = nameLogical;
        fLogicalToMaterialMap[nameLogical] = nameMaterial;
        fLogicalToPhysicalMap[nameLogical].emplace_back(namePhysical);
        fPhysicalToPositionInWorldMap[physicalNewName] = {position.x(), position.y(), position.z()};
        InsertVolumeName(i, physicalNewName);

        if (!fIsAssembly && GetAlternativeNameFromGeant4PhysicalName(namePhysical) != namePhysical) {
            fIsAssembly = true;

            const auto geant4MajorVersionNumber = restG4Metadata->GetGeant4VersionMajor();
            if (geant4MajorVersionNumber < 11) {
                cout << "User geometry consists of assembly which is not supported for this rest / Geant4 "
                        "version combination. Please upgrade to Geant4 11.0.0 or more to use this feature"
                     << endl;
                // exit(1);
            }
        }
    }
}
