
#include "DetectorConstruction.h"

#include <TRestGeant4Geometry.h>
#include <TXMLEngine.h>

#include <G4FieldManager.hh>
#include <G4IonTable.hh>
#include <G4Isotope.hh>
#include <G4MagneticField.hh>
#include <G4Material.hh>
#include <G4SystemOfUnits.hh>
#include <G4UniformMagField.hh>
#include <G4UserLimits.hh>

using namespace std;

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

    G4VPhysicalVolume* worldVolume = parser->GetWorldVolume();

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
            return worldVolume;
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

    auto geometry = restG4Metadata->GetGeometry();
    geometry->LoadGdml(restG4Metadata->Get_GDML_Filename());
    geometry->BuildAssemblyLookupTable(worldVolume);
    geometry->InitializeFromGeant4World(worldVolume);

    for (const auto& activeVolume : geometry->GetActiveVolumes()) {
        cout << " - " << activeVolume << geometry->GetMaxStepSize(activeVolume) << endl;
    }
    for (const auto& activeVolume : geometry->GetActiveVolumes()) {
        G4cout << activeVolume << endl;
        auto physicalVolume = GetPhysicalVolume(activeVolume.Data());

        if (!physicalVolume) {
            G4cout << "DetectorConstruction. Volume " << activeVolume << " is not defined in the geometry"
                   << endl;
            exit(1);
        }

        auto logicalVolume = physicalVolume->GetLogicalVolume();
        if (geometry->GetMaxStepSize(activeVolume) > 0) {
            cout << "step size: " << endl;
            G4cout << "Setting maxStepSize = " << geometry->GetMaxStepSize(activeVolume)
                   << "mm for volume : " << activeVolume << endl;
            logicalVolume->SetUserLimits(new G4UserLimits(geometry->GetMaxStepSize(activeVolume) * mm));
        }
        G4cout << "Activating volume : " << activeVolume << endl;
        restG4Event->AddActiveVolume(activeVolume.Data());
    }

    return worldVolume;
}

G4VPhysicalVolume* DetectorConstruction::GetPhysicalVolume(G4String physVolName) {
    G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
    std::vector<G4VPhysicalVolume*>::const_iterator physVol;
    for (physVol = physVolStore->begin(); physVol != physVolStore->end(); physVol++) {
        if ((*physVol)->GetName() == physVolName) {
            return *physVol;
        }
    }

    return NULL;
}

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
                    G4cout << volumeName << endl;
                }
                attr = xml.GetNextAttr(attr);
            }
        }
        child = xml.GetNext(child);
    }
    return nullptr;
}
TString GetNodeAttribute(TXMLEngine xml, XMLNodePointer_t node, const TString& attributeName) {
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
void AddVolumesRecursively(vector<TString>* container, const vector<TString>& children,
                           map<TString, TString>& nameTable, map<TString, vector<TString>>& childrenTable,
                           const TString& name = "") {
    // G4cout << "called AddVolumesRecursively with name: " << name << endl;
    for (const auto& child : children) {
        // G4cout << "\t" << child << endl;
    }
    if (children.empty()) {
        container->push_back(name);
        // G4cout << "ADDED: " << name << endl;
        return;
    }
    for (const auto& childName : children) {
        AddVolumesRecursively(container, childrenTable[nameTable[childName]], nameTable, childrenTable,
                              name + "_" + childName);
    }
}

void TRestGeant4Geometry::BuildAssemblyLookupTable(const G4VPhysicalVolume* fWorld) {
    // Geometry must be in GDML
    TXMLEngine xml;
    XMLDocPointer_t xmldoc = xml.ParseFile(fGdmlAbsolutePath);
    if (!xmldoc) {
        cout << "TRestGeant4Geometry::BuildAssemblyLookupTable - Failed to open '" << fGdmlAbsolutePath << "'"
             << endl;
        exit(1);
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
    map<TString, TString> lookupTable;
    const int n = int(fWorld->GetLogicalVolume()->GetNoDaughters());
    for (int i = 0; i < n; i++) {
        G4VPhysicalVolume* volume = fWorld->GetLogicalVolume()->GetDaughter(i);
        auto namePhysical = volume->GetName();
        lookupTable[namePhysical] = (*names)[i];
    }
    fGeant4PhysicalToPhysicalMap = lookupTable;
    delete names;
    xml.FreeDoc(xmldoc);
}

void TRestGeant4Geometry::InitializeFromGeant4World(const G4VPhysicalVolume* world) {
    const int n = int(world->GetLogicalVolume()->GetNoDaughters());
    for (int i = 0; i < n; i++) {
        G4VPhysicalVolume* volume = world->GetLogicalVolume()->GetDaughter(i);
        TString namePhysical = (TString)volume->GetName();
        TString physicalLookupAlias = GetPhysicalFromGeant4Physical(TString(namePhysical));
        TString nameLogical = (TString)volume->GetLogicalVolume()->GetName();
        TString nameMaterial = (TString)volume->GetLogicalVolume()->GetMaterial()->GetName();
        auto position = volume->GetTranslation();

        cout << TString::Format(
                    "---> %d - physical: %s (%s)- logical: %s - material: %s - position: (%f, %f, "
                    "%f)",
                    i, namePhysical.Data(), physicalLookupAlias.Data(), nameLogical.Data(),
                    nameMaterial.Data(), position.x(), position.y(), position.z())
             << endl;

        fPhysicalVolumes.emplace_back(namePhysical);
        fPhysicalToLogicalVolumeMap[namePhysical] = nameLogical;
        fLogicalToMaterialMap[nameLogical] = nameMaterial;
        fPhysicalToPositionInWorldMap[namePhysical] = {position.x(), position.y(), position.z()};
    }
    if (fAllVolumesActive || fActiveVolumes.empty()) {
        for (const auto& physicalVolume : fPhysicalVolumes) {
            InsertActiveVolume(physicalVolume, fDefaultStorageChance, fDefaultMaxStepSize);
        }
    }
}