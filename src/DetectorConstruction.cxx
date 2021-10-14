
#include "DetectorConstruction.h"

#include <TXMLEngine.h>
#include <spdlog/spdlog.h>

#include <G4FieldManager.hh>
#include <G4IonTable.hh>
#include <G4Isotope.hh>
#include <G4MagneticField.hh>
#include <G4Material.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UniformMagField.hh>
#include <G4UserLimits.hh>

#include "GlobalManager.h"
#include "SensitiveDetector.h"

using namespace std;

DetectorConstruction::DetectorConstruction()
    : fRestGeant4Metadata(GlobalManager::Instance()->GetRestGeant4Metadata()) {
    G4cout << "Detector Construction" << G4endl;
    parser = new G4GDMLParser();
    fGeometryFilename = fRestGeant4Metadata->Get_GDML_Filename();
}

DetectorConstruction::~DetectorConstruction() { delete parser; }

G4VPhysicalVolume* DetectorConstruction::Construct() {
    cout << "Isotope table " << endl;
    cout << *(G4Isotope::GetIsotopeTable()) << endl;

    G4cout << "Producing geometry" << G4endl;

    // Reading the geometry
    TString geometryFile = fRestGeant4Metadata->Get_GDML_Filename();

    char originDirectory[256];
    sprintf(originDirectory, "%s", getenv("PWD"));
    auto separatePathAndName =
        TRestTools::SeparatePathAndName((string)fRestGeant4Metadata->Get_GDML_Filename());
    chdir(separatePathAndName.first.c_str());

    parser->Read(separatePathAndName.second, false);

    chdir(originDirectory);

    fWorld = parser->GetWorldVolume();

    // TODO : Take the name of the sensitive volume and use it here to define its
    // StepSize
    string SensVol = (string)fRestGeant4Metadata->GetSensitiveVolume();
    G4VPhysicalVolume* _vol = GetPhysicalVolume(SensVol);
    if (!_vol) {
        G4cout << "RESTG4 error. Sensitive volume " << SensVol << " does not exist in geometry!!" << G4endl;
        PrintGeometryInfo();
        G4cout << "RESTG4 error. Please, review geometry! Presh a key to crash!!" << G4endl;
        getchar();
        // We need to produce a clean exit at this point
        exit(1);
    }

    Double_t mx = fRestGeant4Metadata->GetMagneticField().X() * tesla;
    Double_t my = fRestGeant4Metadata->GetMagneticField().Y() * tesla;
    Double_t mz = fRestGeant4Metadata->GetMagneticField().Z() * tesla;

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
        G4cout << "Sensitive volume properties" << G4endl;
        G4cout << "==============" << G4endl;
        G4cout << "Sensitive volume name: " << mat->GetName() << G4endl;
        G4cout << "Sensitive volume temperature: " << mat->GetTemperature() << G4endl;
        G4cout << "Sensitive volume density: " << mat->GetDensity() / (g / cm3) << " g/cm3" << G4endl;
    } else {
        cout << "ERROR : Logical volume for sensitive \"" << SensVol << "\" not found!" << endl;
    }

    // Getting generation volume
    string GenVol = (string)fRestGeant4Metadata->GetGeneratedFrom();
    cout << "Generated from volume: " << GenVol << endl;
    string type = (string)fRestGeant4Metadata->GetGeneratorType();
    cout << "Generator type: " << type << endl;

    // TODO if we do not find the volume given in the config inside the geometry
    // we should RETURN error
    if (type == "volume" && GenVol != "Not defined") {
        G4VPhysicalVolume* pVol = GetPhysicalVolume(GenVol);
        if (pVol == nullptr) {
            cout << "ERROR : The generator volume was not found in the geometry" << endl;
            exit(1);
            return fWorld;
        }

        fGeneratorTranslation = pVol->GetTranslation();

        // We set in TRestGeant4Metadata the center of the generator. If it is a point
        // we just want the value from the config file.
        // TODO : make this kind of keyword comparisons case insensitive?
        if (type == "surface" || type == "volume") {
            fRestGeant4Metadata->SetGeneratorPosition(fGeneratorTranslation.x(), fGeneratorTranslation.y(),
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

    for (int id = 0; id < fRestGeant4Metadata->GetNumberOfActiveVolumes(); id++) {
        TString actVolName = fRestGeant4Metadata->GetActiveVolumeName(id);
        G4VPhysicalVolume* pVol = GetPhysicalVolume((G4String)actVolName);

        cout << "Activating volume : " << actVolName << endl;
        // restG4Event->AddActiveVolume((string)actVolName);

        if (!pVol) {
            cout << "DetectorConstruction. Volume " << actVolName << " is not defined in the geometry"
                 << endl;
            // exit(1); TODO FIX THIS FOR ASSEMBLY
        } else {
            G4LogicalVolume* lVol = pVol->GetLogicalVolume();
            if (fRestGeant4Metadata->GetMaxStepSize(actVolName) > 0) {
                G4cout << "Setting maxStepSize = " << fRestGeant4Metadata->GetMaxStepSize(actVolName)
                       << "mm for volume : " << actVolName << G4endl;
                lVol->SetUserLimits(new G4UserLimits(fRestGeant4Metadata->GetMaxStepSize(actVolName) * mm));
            }
        }
    }

    cout << "Detector constructed : " << fWorld << endl;

    return fWorld;
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

void DetectorConstruction::ConstructSDandField() {
    auto restGeant4Metadata = GlobalManager::Instance()->GetRestGeant4Metadata();
    TString sensitiveVolumeName = restGeant4Metadata->GetSensitiveVolume();

    if (sensitiveVolumeName.IsNull()) {
        return;
    }

    G4SDManager* SDManager = G4SDManager::GetSDMpointer();

    /*
    if (!fSensitiveLogicalVolumeNames.empty()) {
        for (const auto& fSensitiveLogicalVolumeName : fSensitiveLogicalVolumeNames) {
            G4LogicalVolume* logicalVolume =
                G4LogicalVolumeStore::GetInstance()->GetVolume(fSensitiveLogicalVolumeName);
            if (!logicalVolume) {
                // PrintGeometryInfo();
                cout
                    << "Trying to attach a sensitive detector to logical volume '{}', but this volume is not "
                       "found in store."
                    << fSensitiveLogicalVolumeName << endl;
                exit(1);
            }
            G4VSensitiveDetector* sensitiveDetector = new SensitiveDetector(fSensitiveLogicalVolumeName);
            SDManager->AddNewDetector(sensitiveDetector);
            logicalVolume->SetSensitiveDetector(sensitiveDetector);
            spdlog::info("Attaching Sensitive detector '{}' to logical volume '{}'",
                         sensitiveDetector->GetName(), fSensitiveLogicalVolumeName);
            for (int i = 0; i < int(fWorld->GetLogicalVolume()->GetNoDaughters()); i++) {
                auto physicalVolume = fWorld->GetLogicalVolume()->GetDaughter(i);
                const string name = physicalVolume->GetLogicalVolume()->GetName();
                if (name == fSensitiveLogicalVolumeName)
                    spdlog::info("---> Physical volume: {}", physicalVolume->GetName());
            }
        }
    }
    */

    vector<string> fSensitivePhysicalVolumeNames = {(string)(sensitiveVolumeName)};
    if (!fSensitivePhysicalVolumeNames.empty()) {
        for (const auto& fSensitivePhysicalVolumeName : fSensitivePhysicalVolumeNames) {
            G4VPhysicalVolume* physicalVolume =
                G4PhysicalVolumeStore::GetInstance()->GetVolume(fSensitivePhysicalVolumeName);
            if (!physicalVolume) {
                // PrintGeometryInfo();
                cout << "Trying to attach a sensitive detector to the logical volume of physical volume '{}'"
                        ", but this physical volume is not found in store."
                     << fSensitivePhysicalVolumeName << endl;
                exit(1);
            }
            G4LogicalVolume* logicalVolume = physicalVolume->GetLogicalVolume();
            G4VSensitiveDetector* sensitiveDetector = new SensitiveDetector(logicalVolume->GetName());
            SDManager->AddNewDetector(sensitiveDetector);
            logicalVolume->SetSensitiveDetector(sensitiveDetector);
        }
    }
}

void DetectorConstruction::PrintGeometryInfo() {
    spdlog::info("DetectorConstruction::PrintGeometryInfo - Begin");
    const int n = int(fWorld->GetLogicalVolume()->GetNoDaughters());
    for (int i = 0; i < n; i++) {
        G4VPhysicalVolume* volume = fWorld->GetLogicalVolume()->GetDaughter(i);
        auto namePhysical = volume->GetName();
        auto nameLogical = volume->GetLogicalVolume()->GetName();
        auto nameMaterial = volume->GetLogicalVolume()->GetMaterial()->GetName();
        auto position = volume->GetTranslation();
        auto physicalLookupAlias = GlobalManager::Instance()->GetVolumeFromLookupTable(namePhysical);
        spdlog::info(
            "---> {} - physical: {} ({})- logical: {} - material: {} - position: ({:03.2f}, {:03.2f}, "
            "{:03.2f})",
            i, namePhysical, physicalLookupAlias, nameLogical, nameMaterial, position.x(), position.y(),
            position.z());
    }
    spdlog::info("DetectorConstruction::PrintGeometryInfo - End");
}

/*
XMLNodePointer_t FindChildByName(TXMLEngine xml, XMLNodePointer_t node, const TString& name) {
    XMLNodePointer_t child = xml.GetChild(node);
    while (child) {
        TString childName = xml.GetNodeName(child);
        if (childName.EqualTo(name)) {
            return child;
        }
        child = xml.GetNext(child);
    }
    return nullptr;
}

XMLNodePointer_t FindVolumeOrAssemblyByName(TXMLEngine xml, XMLNodePointer_t node, const TString& name) {
    XMLNodePointer_t child = xml.GetChild(node);
    while (child) {
        TString childName = xml.GetNodeName(child);
        if (childName.EqualTo("volume") || childName.EqualTo("assembly")) {
            XMLAttrPointer_t attr = xml.GetFirstAttr(child);
            while (attr) {
                if (TString(xml.GetAttrName(attr)).EqualTo("name")) {
                    TString volumeName = xml.GetAttrValue(attr);
                    cout << volumeName << endl;
                }
                attr = xml.GetNextAttr(attr);
            }
        }
        child = xml.GetNext(child);
    }

    return nullptr;
}

TString GetNodeAttribute(TXMLEngine xml, XMLNodePointer_t node, TString attributeName) {
    XMLAttrPointer_t attr = xml.GetFirstAttr(node);
    while (attr) {
        if (TString(xml.GetAttrName(attr)).EqualTo(attributeName)) {
            TString refName = xml.GetAttrValue(attr);
            return refName;
        }
        attr = xml.GetNextAttr(attr);
    }
    return TString();
}

void AddVolumesRecursively(vector<TString>* container, vector<TString> children,
                           map<TString, TString>& nameTable, map<TString, vector<TString>>& childrenTable,
                           const TString& name = "") {
    cout << "called AddVolumesRecursively with name: " << name << endl;
    for (const auto& child : children) {
        cout << "\t" << child << endl;
    }

    if (children.size() == 0) {
        container->push_back(name);
        cout << "ADDED: " << name << endl;
        return;
    }
    for (const auto& childName : children) {
        AddVolumesRecursively(container, childrenTable[nameTable[childName]], nameTable, childrenTable,
                              name + "_" + childName);
    }
}

void DetectorConstruction::BuildAssemblyLookupTable() {
    // spdlog::info("DetectorConstruction::BuildAssemblyLookupTable - Begin");
    // Only if geometry is in gdml

    if (fGeometryFilename.substr(fGeometryFilename.length() - 5) != ".gdml") {
        // spdlog::debug("DetectorConstruction::BuildAssemblyLookupTable - Geometry not GDML, doing nothing");
        return;
    }

    TXMLEngine xml;
    XMLDocPointer_t xmldoc = xml.ParseFile(fGeometryFilename.c_str());
    if (!xmldoc) {
        // spdlog::warn("DetectorConstruction::BuildAssemblyLookupTable - Problems reading GDML file '{}'",
        //              fGeometryFilename);
        return;
    }

    map<TString, TString> nameTable;
    map<TString, vector<TString>> childrenTable;

    XMLNodePointer_t mainNode = xml.DocGetRootElement(xmldoc);
    XMLNodePointer_t structure = FindChildByName(xml, mainNode, "structure");
    XMLNodePointer_t world = FindVolumeOrAssemblyByName(xml, structure, "world");

    XMLNodePointer_t child = xml.GetChild(structure);
    while (child) {
        TString name = xml.GetNodeName(child);
        TString volumeName = GetNodeAttribute(xml, child, "name");
        auto physicalVolumeNode = xml.GetChild(child);
        childrenTable[volumeName] = {};

        while (physicalVolumeNode) {
            auto physicalVolumeName = GetNodeAttribute(xml, physicalVolumeNode, "name");
            auto volumeRefNode = xml.GetChild(physicalVolumeNode);
            while (volumeRefNode) {
                TString volumeRefNodeName = xml.GetNodeName(volumeRefNode);
                if (volumeRefNodeName.EqualTo("volumeref")) {
                    TString refName = GetNodeAttribute(xml, volumeRefNode, "ref");
                    nameTable[physicalVolumeName] = refName;
                    childrenTable[volumeName].push_back(physicalVolumeName);
                }
                volumeRefNode = xml.GetNext(volumeRefNode);
            }
            physicalVolumeNode = xml.GetNext(physicalVolumeNode);
        }

        child = xml.GetNext(child);
    }

    auto names = new vector<TString>();
    for (const auto& topName : childrenTable["world"]) {
        auto children = childrenTable[nameTable[topName]];
        AddVolumesRecursively(names, children, nameTable, childrenTable, topName);
    }

    map<string, string> lookupTable;

    const int n = int(fWorld->GetLogicalVolume()->GetNoDaughters());
    for (int i = 0; i < n; i++) {
        G4VPhysicalVolume* volume = fWorld->GetLogicalVolume()->GetDaughter(i);
        auto namePhysical = volume->GetName();
        lookupTable[namePhysical] = (*names)[i];
    }

    GlobalManager::Instance()->SetVolumeLookupTable(lookupTable);

    delete names;
    xml.FreeDoc(xmldoc);

    // spdlog::info("DetectorConstruction::BuildAssemblyLookupTable - End");
}
*/