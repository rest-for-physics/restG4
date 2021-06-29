
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <TRestGeant4PhysicsLists.h>

#include <G4EmConfigurator.hh>
#include <G4VModularPhysicsList.hh>
#include <globals.hh>

using namespace std;

class G4VPhysicsConstructor;

class PhysicsList : public G4VModularPhysicsList {
   public:
    PhysicsList();
    PhysicsList(TRestGeant4PhysicsLists* restPhysicsLists);
    ~PhysicsList();

   protected:
    // Construct particle and physics
    virtual void InitializePhysicsLists();
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    virtual void SetCuts();

   private:
    G4EmConfigurator em_config;

    G4VPhysicsConstructor* fEmPhysicsList;
    G4VPhysicsConstructor* fDecPhysicsList;
    G4VPhysicsConstructor* fRadDecPhysicsList;
    std::vector<G4VPhysicsConstructor*> fHadronPhys;

    TRestGeant4PhysicsLists* restPhysList;
};

#endif
