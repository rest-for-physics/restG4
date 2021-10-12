//
// Created by lobis on 12/10/2021.
//

#ifndef REST_ACTIONINITIALIZATION_H
#define REST_ACTIONINITIALIZATION_H

#include <G4VUserActionInitialization.hh>

class DetectorConstruction;
class G4VSteppingVerbose;

class ActionInitialization : public G4VUserActionInitialization {
   public:
    ActionInitialization();
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

    virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;
};

#endif  // REST_ACTIONINITIALIZATION_H
