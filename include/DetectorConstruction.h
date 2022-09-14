
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>

#include <G4GDMLParser.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VUserDetectorConstruction.hh>

class SimulationManager;

class DetectorConstruction : public G4VUserDetectorConstruction {
   private:
    SimulationManager* fSimulationManager;

    G4GDMLParser* fGdmlParser;
    G4VSolid* fGeneratorSolid;
    G4ThreeVector fGeneratorTranslation;

    Double_t fBoundBoxXMin, fBoundBoxXMax, fBoundBoxYMin, fBoundBoxYMax, fBoundBoxZMin, fBoundBoxZMax;

   public:
    G4VPhysicalVolume* GetPhysicalVolume(const G4String& physVolName) const;
    inline G4VSolid* GetGeneratorSolid() const { return fGeneratorSolid; }
    inline G4ThreeVector GetGeneratorTranslation() const { return fGeneratorTranslation; }

    inline Double_t GetBoundBoxXMin() const { return fBoundBoxXMin; }
    inline Double_t GetBoundBoxXMax() const { return fBoundBoxXMax; }
    inline Double_t GetBoundBoxYMin() const { return fBoundBoxYMin; }
    inline Double_t GetBoundBoxYMax() const { return fBoundBoxYMax; }
    inline Double_t GetBoundBoxZMin() const { return fBoundBoxZMin; }
    inline Double_t GetBoundBoxZMax() const { return fBoundBoxZMax; }

    DetectorConstruction(SimulationManager*);
    ~DetectorConstruction();

   public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    friend class TRestGeant4GeometryInfo;
};
#endif
