
#ifndef REST_SENSITIVEDETECTOR_H
#define REST_SENSITIVEDETECTOR_H

#include <G4VSensitiveDetector.hh>

class G4TouchableHistory;
class SimulationManager;

class SensitiveDetector : public G4VSensitiveDetector {
   public:
    explicit SensitiveDetector(SimulationManager* simulationManager, const G4String& name);
    virtual ~SensitiveDetector() = default;

    virtual void Initialize(G4HCofThisEvent*) {}
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);

   private:
    SimulationManager* fSimulationManager;
};
#endif  // REST_SENSITIVEDETECTOR_H
