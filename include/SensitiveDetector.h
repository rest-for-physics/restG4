//
// Created by lobis on 12/10/2021.
//

#ifndef REST_SENSITIVEDETECTOR_H
#define REST_SENSITIVEDETECTOR_H

#include <G4VSensitiveDetector.hh>
#include <string>

class OutputManager;

class G4TouchableHistory;

class SensitiveDetector : public G4VSensitiveDetector {
   public:
    explicit SensitiveDetector(const G4String& name);
    virtual ~SensitiveDetector() = default;

    virtual void Initialize(G4HCofThisEvent*) {}
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);

   private:
    OutputManager* fOutputManager;
};

#endif  // REST_SENSITIVEDETECTOR_H
