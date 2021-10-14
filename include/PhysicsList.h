
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <TRestGeant4PhysicsLists.h>

#include <G4EmConfigurator.hh>
#include <G4VModularPhysicsList.hh>
#include <globals.hh>

class G4VPhysicsConstructor;

class PhysicsList : public G4VModularPhysicsList {
   public:
    PhysicsList(G4int verbosity = 0);
    explicit PhysicsList(TRestGeant4PhysicsLists* restPhysicsLists);
    ~PhysicsList() override;

   protected:
    // Construct particle and physics
    virtual void InitializePhysicsLists();
    void ConstructParticle() override;
    void ConstructProcess() override;
    void SetCuts() override;
    void SetCutValue(G4double, const G4String&);

   private:
    G4EmConfigurator emConfigurator;

    G4VPhysicsConstructor* fEmPhysicsList;
    G4VPhysicsConstructor* fDecayPhysicsList;
    G4VPhysicsConstructor* fRadioactiveDecayPhysicsList;
    std::vector<G4VPhysicsConstructor*> fHadronPhys;

    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists;
};

#endif
