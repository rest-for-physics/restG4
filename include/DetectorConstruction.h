
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>

#include <G4GDMLParser.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VUserDetectorConstruction.hh>

class DetectorConstruction : public G4VUserDetectorConstruction {
   private:
    G4VPhysicalVolume* fWorld;

    TRestGeant4Metadata* fRestGeant4Metadata;

    G4GDMLParser* parser;
    G4VSolid* generatorSolid;
    G4ThreeVector fGeneratorTranslation;

    string fGeometryFilename;

    Double_t boundBox_xMin, boundBox_xMax;
    Double_t boundBox_yMin, boundBox_yMax;
    Double_t boundBox_zMin, boundBox_zMax;

   public:
    void ConstructSDandField() override;
    void BuildAssemblyLookupTable();

    G4GDMLParser* GetGeometry() { return parser; }
    G4VPhysicalVolume* GetPhysicalVolume(G4String physVolName);
    G4VSolid* GetGeneratorSolid() { return generatorSolid; }
    G4ThreeVector GetGeneratorTranslation() { return fGeneratorTranslation; }

    Double_t GetBoundingX_min() { return boundBox_xMin; }
    Double_t GetBoundingX_max() { return boundBox_xMax; }

    Double_t GetBoundingY_min() { return boundBox_yMin; }
    Double_t GetBoundingY_max() { return boundBox_yMax; }

    Double_t GetBoundingZ_min() { return boundBox_zMin; }
    Double_t GetBoundingZ_max() { return boundBox_zMax; }

    DetectorConstruction();
    ~DetectorConstruction();

   public:
    G4VPhysicalVolume* Construct();
    void PrintGeometryInfo();
};
#endif
