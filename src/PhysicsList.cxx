
#include "PhysicsList.h"

#include <CLHEP/Units/PhysicalConstants.h>

#include <G4BetheBlochIonGasModel.hh>
#include <G4BraggIonGasModel.hh>
#include <G4DecayPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4EmLivermorePhysics.hh>
#include <G4EmParameters.hh>
#include <G4EmPenelopePhysics.hh>
#include <G4EmStandardPhysics_option3.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4HadronElasticPhysicsHP.hh>
#include <G4HadronPhysicsQGSP_BIC_HP.hh>
#include <G4IonBinaryCascadePhysics.hh>
#include <G4IonFluctuations.hh>
#include <G4IonParametrisedLossModel.hh>
#include <G4IonTable.hh>
#include <G4LossTableManager.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleTypes.hh>
#include <G4ProcessManager.hh>
#include <G4ProductionCuts.hh>
#include <G4RadioactiveDecay.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4StepLimiter.hh>
#include <G4StoppingPhysics.hh>
#include <G4SystemOfUnits.hh>
#include <G4UAtomicDeexcitation.hh>
#include <G4UnitsTable.hh>
#include <G4UniversalFluctuation.hh>

#include "Particles.h"

Int_t emCounter = 0;

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
    cout << "restG4. PhysicsList. Wrong constructor!!" << endl;
}

PhysicsList::PhysicsList(TRestGeant4PhysicsLists* physicsLists) : G4VModularPhysicsList() {
    // add new units for radioActive decays
    const G4double minute = 60 * second;
    const G4double hour = 60 * minute;
    const G4double day = 24 * hour;
    const G4double year = 365 * day;
    new G4UnitDefinition("minute", "min", "Time", minute);
    new G4UnitDefinition("hour", "h", "Time", hour);
    new G4UnitDefinition("day", "d", "Time", day);
    new G4UnitDefinition("year", "y", "Time", year);

    defaultCutValue = 0.1 * mm;

    restPhysList = physicsLists;
    G4LossTableManager::Instance();
    // fix lower limit for cut
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(
        restPhysList->GetMinimumEnergyProductionCuts() * keV,
        restPhysList->GetMaximumEnergyProductionCuts() * keV);
    defaultCutValue = 0.1 * mm;

    fEmPhysicsList = nullptr;
    fDecPhysicsList = nullptr;
    fRadDecPhysicsList = nullptr;

    InitializePhysicsLists();
}

PhysicsList::~PhysicsList() {
    delete fEmPhysicsList;

    delete fDecPhysicsList;
    delete fRadDecPhysicsList;

    for (auto hadronicPhysicsList : fHadronPhys) {
        delete hadronicPhysicsList;
    }
}

void PhysicsList::InitializePhysicsLists() {
    // Decay physics and all particles
    if (restPhysList->FindPhysicsList("G4DecayPhysics") >= 0)
        fDecPhysicsList = new G4DecayPhysics();
    else if (restPhysList->GetVerboseLevel() >= REST_Debug) {
        G4cout << "restG4. PhysicsList. G4DecayPhysics is not enabled!!" << G4endl;
    }

    // RadioactiveDecay physicsList
    if (restPhysList->FindPhysicsList("G4RadioactiveDecayPhysics") >= 0)
        fRadDecPhysicsList = new G4RadioactiveDecayPhysics();
    else if (restPhysList->GetVerboseLevel() >= REST_Debug) {
        G4cout << "restG4. PhysicsList. G4RadioactiveDecayPhysics is not enabled!!" << G4endl;
    }

    // Electromagnetic physicsList
    if (restPhysList->FindPhysicsList("G4EmLivermorePhysics") >= 0) {
        if (fEmPhysicsList == nullptr) fEmPhysicsList = new G4EmLivermorePhysics();
        emCounter++;
    }

    if (restPhysList->FindPhysicsList("G4EmPenelopePhysics") >= 0) {
        if (fEmPhysicsList == nullptr) fEmPhysicsList = new G4EmPenelopePhysics();
        emCounter++;
    }

    if (restPhysList->FindPhysicsList("G4EmStandardPhysics_option3") >= 0) {
        if (fEmPhysicsList == nullptr) fEmPhysicsList = new G4EmStandardPhysics_option3();
        emCounter++;
    }

    if (restPhysList->FindPhysicsList("G4EmStandardPhysics_option4") >= 0) {
        if (fEmPhysicsList == nullptr) fEmPhysicsList = new G4EmStandardPhysics_option4();
        emCounter++;
    }

    if (restPhysList->GetVerboseLevel() >= REST_Essential && emCounter == 0) {
        G4cout << "REST WARNING : No electromagenetic physics list has been enabled!!" << G4endl;
    }

    if (emCounter > 1) {
        G4cout << "REST ERROR: restG4. PhysicsList. More than 1 EM PhysicsList "
                  "enabled."
               << G4endl;
        exit(1);
    }

    // Hadronic PhysicsList
    if (restPhysList->FindPhysicsList("G4HadronPhysicsQGSP_BIC_HP") >= 0)
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BIC_HP());

    if (restPhysList->FindPhysicsList("G4IonBinaryCascadePhysics") >= 0)
        fHadronPhys.push_back(new G4IonBinaryCascadePhysics());

    if (restPhysList->FindPhysicsList("G4HadronElasticPhysicsHP") >= 0)
        fHadronPhys.push_back(new G4HadronElasticPhysicsHP());

    if (restPhysList->FindPhysicsList("G4NeutronTrackingCut") >= 0)
        fHadronPhys.push_back(new G4NeutronTrackingCut());

    if (restPhysList->FindPhysicsList("G4EmExtraPhysics") >= 0) fHadronPhys.push_back(new G4EmExtraPhysics());

    G4cout << "Number of hadronic physics lists added " << fHadronPhys.size() << G4endl;
}

void PhysicsList::ConstructParticle() {
    // pseudo-particles
    G4Geantino::GeantinoDefinition();

    // particles defined in PhysicsLists
    if (fDecPhysicsList) fDecPhysicsList->ConstructParticle();

    if (fEmPhysicsList) fEmPhysicsList->ConstructParticle();

    if (fRadDecPhysicsList) fRadDecPhysicsList->ConstructParticle();

    for (size_t i = 0; i < fHadronPhys.size(); i++) fHadronPhys[i]->ConstructParticle();
}

void PhysicsList::ConstructProcess() {
    AddTransportation();

    // Electromagnetic physics list
    if (fEmPhysicsList) {
        fEmPhysicsList->ConstructProcess();
        em_config.AddModels();
    }

    // Decay physics list
    if (fDecPhysicsList) fDecPhysicsList->ConstructProcess();

    // Radioactive decay
    if (fRadDecPhysicsList) fRadDecPhysicsList->ConstructProcess();

    // Hadronic physics lists
    for (size_t i = 0; i < fHadronPhys.size(); i++) fHadronPhys[i]->ConstructProcess();

    if (restPhysList->FindPhysicsList("G4RadioactiveDecay")) {
        auto radioactiveDecay = new G4RadioactiveDecay();

        radioactiveDecay->SetThresholdForVeryLongDecayTime(nanosecond);

        // Setting Internal Conversion (ICM) option.
        if (restPhysList->GetPhysicsListOptionValue("G4RadioactiveDecay", "ICM") == "true")
            radioactiveDecay->SetICM(true);  // Internal Conversion
        else if (restPhysList->GetPhysicsListOptionValue("G4RadioactiveDecay", "ICM") == "false")
            radioactiveDecay->SetICM(false);  // Internal Conversion
        else if (restPhysList->GetVerboseLevel() >= REST_Essential)
            G4cout << "REST WARNING. restG4. PhysicsList. G4RadioactiveDecay. Option "
                      "ICM not defined."
                   << G4endl;

        // Enabling electron re-arrangement (ARM) option.
        if (restPhysList->GetPhysicsListOptionValue("G4RadioactiveDecay", "ARM") == "true")
            radioactiveDecay->SetARM(true);  // Internal Conversion
        else if (restPhysList->GetPhysicsListOptionValue("G4RadioactiveDecay", "ARM") == "false")
            radioactiveDecay->SetARM(false);  // Internal Conversion
        else if (restPhysList->GetVerboseLevel() >= REST_Essential)
            G4cout << "REST WARNING. restG4. PhysicsList. G4RadioactiveDecay. Option "
                      "ARM not defined."
                   << G4endl;
    }

/*
   G4ScreenedNuclearRecoil* nucr = new G4ScreenedNuclearRecoil();
   G4double energyLimit = 100.*MeV;
   nucr->SetMaxEnergyForScattering(energyLimit);
   ph->RegisterProcess( nucr, G4GenericIon::GenericIon());
   */

/* theParticleIterator->reset();
while ((*theParticleIterator)())
{
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String partname = particle->GetParticleName();
    if(partname == "alpha" || partname == "He3" || partname == "GenericIon") {
        G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
        G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
        G4double eth = 2.*MeV*particle->GetPDGMass()/CLHEP::proton_mass_c2;
        em_config.SetExtraEmModel(partname,"braggIoni",mod1,"",0.0,eth,
                new G4IonFluctuations());
        em_config.SetExtraEmModel(partname,"betheIoni",mod2,"",eth,100*TeV,
                new G4UniversalFluctuation());

    }
} */

// Requires Geant4 version higher than 10.2.9. Defined at CMakeLists.
#ifdef G4104
    auto theParticleIterator = GetParticleIterator();

    // to implement UserLimits to StepSize inside the gas
    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
        G4ParticleDefinition* particle = theParticleIterator->value();
        const auto& particleName = particle->GetParticleName();
        G4ProcessManager* processManager = particle->GetProcessManager();

        if (particleName == "e-")
            processManager->AddDiscreteProcess(new G4StepLimiter("e-Step"));
        else if (particleName == "e+")
            processManager->AddDiscreteProcess(new G4StepLimiter("e+Step"));

        if (particleName == "mu-")
            processManager->AddDiscreteProcess(new G4StepLimiter("mu-Step"));
        else if (particleName == "mu+")
            processManager->AddDiscreteProcess(new G4StepLimiter("mu+Step"));
    }

    // There might be a better way to do this
    for (int Z = 1; Z <= 40; Z++)
        for (int A = 2 * Z; A <= 3 * Z; A++) {
            for (unsigned int n = 0; n < restPhysList->GetIonStepList().size(); n++) {
                if (restPhysList->GetIonStepList()[n] == G4IonTable::GetIonTable()->GetIonName(Z, A)) {
                    G4ParticleDefinition* particle = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
                    G4String particle_name = G4IonTable::GetIonTable()->GetIonName(Z, A, 0);
                    cout << "Found ion: " << particle_name << " Z " << Z << " A " << A << endl;
                    G4ProcessManager* processManager = particle->GetProcessManager();
                    processManager->AddDiscreteProcess(new G4StepLimiter("ionStep"));
                }
            }
        }
#endif  // G4104
}

void PhysicsList::SetCuts() {
    SetCutsWithDefault();

    SetCutValue(restPhysList->GetCutForGamma() * mm, "gamma");
    SetCutValue(restPhysList->GetCutForElectron() * mm, "e-");
    SetCutValue(restPhysList->GetCutForPositron() * mm, "e+");
    SetCutValue(restPhysList->GetCutForMuon() * mm, "mu+");
    SetCutValue(restPhysList->GetCutForMuon() * mm, "mu-");
    SetCutValue(restPhysList->GetCutForNeutron() * mm, "neutron");
}
