
#include "PhysicsList.h"

#include <CLHEP/Units/PhysicalConstants.h>
#include <spdlog/spdlog.h>

#include <G4BetheBlochIonGasModel.hh>
#include <G4BraggIonGasModel.hh>
#include <G4DecayPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4EmLivermorePhysics.hh>
#include <G4EmParameters.hh>
#include <G4EmPenelopePhysics.hh>
#include <G4EmProcessOptions.hh>
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

using namespace std;

Int_t emCounter = 0;

PhysicsList::PhysicsList(G4int verbosity) : G4VModularPhysicsList() {
    spdlog::debug("PhysicsList::PhysicsList");
    SetVerboseLevel(verbosity);
    // add new units for radioActive decays
    //
    const G4double minute = 60 * second;
    const G4double hour = 60 * minute;
    const G4double day = 24 * hour;
    const G4double year = 365 * day;
    new G4UnitDefinition("minute", "min", "Time", minute);
    new G4UnitDefinition("hour", "h", "Time", hour);
    new G4UnitDefinition("day", "d", "Time", day);
    new G4UnitDefinition("year", "y", "Time", year);
}

PhysicsList::PhysicsList(TRestGeant4PhysicsLists* physicsLists) : PhysicsList() {
    fRestGeant4PhysicsLists = physicsLists;
    G4LossTableManager::Instance();

    fEmPhysicsList = nullptr;
    fDecayPhysicsList = nullptr;
    fRadioactiveDecayPhysicsList = nullptr;

    InitializePhysicsLists();
}

PhysicsList::~PhysicsList() { return; }

void PhysicsList::InitializePhysicsLists() {
    // Decay physics and all particles
    if (fRestGeant4PhysicsLists->FindPhysicsList("G4DecayPhysics") >= 0)
        fDecayPhysicsList = new G4DecayPhysics();
    else if (fRestGeant4PhysicsLists->GetVerboseLevel() >= REST_Debug) {
        G4cout << "restG4. PhysicsList. G4DecayPhysics is not enabled!!" << G4endl;
    }

    // RadioactiveDecay physicsList
    if (fRestGeant4PhysicsLists->FindPhysicsList("G4RadioactiveDecayPhysics") >= 0)
        fRadioactiveDecayPhysicsList = new G4RadioactiveDecayPhysics();
    else if (fRestGeant4PhysicsLists->GetVerboseLevel() >= REST_Debug) {
        G4cout << "restG4. PhysicsList. G4RadioactiveDecayPhysics is not enabled!!" << G4endl;
    }

    // Electromagnetic physicsList
    if (fRestGeant4PhysicsLists->FindPhysicsList("G4EmLivermorePhysics") >= 0) {
        if (!fEmPhysicsList) fEmPhysicsList = new G4EmLivermorePhysics(verboseLevel);
        emCounter++;
    }

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4EmPenelopePhysics") >= 0) {
        if (!fEmPhysicsList) fEmPhysicsList = new G4EmPenelopePhysics(verboseLevel);
        emCounter++;
    }

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4EmStandardPhysics_option3") >= 0) {
        if (!fEmPhysicsList) fEmPhysicsList = new G4EmStandardPhysics_option3(verboseLevel);
        emCounter++;
    }

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4EmStandardPhysics_option4") >= 0) {
        if (!fEmPhysicsList) fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);
        emCounter++;
    }

    if (emCounter != 1) {
        if (emCounter == 0) {
            spdlog::warn(
                "PhysicsList::InitializePhysicsLists ---> no electromagnetic physics list has been "
                "selected!");
        } else {
            spdlog::warn(
                "PhysicsList::InitializePhysicsLists ---> more than one electromagnetic physics list has "
                "been selected, only one can be selected");
        }
        exit(1);
    }

    // Hadronic PhysicsList
    if (fRestGeant4PhysicsLists->FindPhysicsList("G4HadronPhysicsQGSP_BIC_HP") >= 0) {
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BIC_HP());
    }

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4IonBinaryCascadePhysics") >= 0) {
        fHadronPhys.push_back(new G4IonBinaryCascadePhysics());
    }

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4HadronElasticPhysicsHP") >= 0) {
        fHadronPhys.push_back(new G4HadronElasticPhysicsHP());
    }

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4NeutronTrackingCut") >= 0) {
        fHadronPhys.push_back(new G4NeutronTrackingCut());
    }

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4EmExtraPhysics") >= 0) {
        fHadronPhys.push_back(new G4EmExtraPhysics());
    }

    spdlog::info("PhysicsList::InitializePhysicsLists ---> Number of hadronic physics lists added: {}",
                 fHadronPhys.size());

    RegisterPhysics(fEmPhysicsList);
    RegisterPhysics(fDecayPhysicsList);
    RegisterPhysics(fRadioactiveDecayPhysicsList);

    for (auto& fHadronPhy : fHadronPhys) {
        RegisterPhysics(fHadronPhy);
    }
}

void PhysicsList::ConstructParticle() { G4VModularPhysicsList::ConstructParticle(); }

void PhysicsList::ConstructProcess() {
    G4VModularPhysicsList::ConstructProcess();

    // TODO: we need to do some refactoring here, most of these lines are redundant (maybe all?)
    AddTransportation();
    // Electromagnetic physics list
    if (fEmPhysicsList) {
        fEmPhysicsList->ConstructProcess();
        emConfigurator.AddModels();
        G4EmProcessOptions emOptions;
        emOptions.SetFluo(true);   // To activate deexcitation processes and fluorescence
        emOptions.SetAuger(true);  // To activate Auger effect if deexcitation is activated
        emOptions.SetPIXE(true);   // To activate Particle Induced X-Ray Emission (PIXE)
    }

    // Decay physics list
    if (fDecayPhysicsList) fDecayPhysicsList->ConstructProcess();

    // Radioactive decay
    if (fRadioactiveDecayPhysicsList) fRadioactiveDecayPhysicsList->ConstructProcess();

    // hadronic physics lists
    for (size_t i = 0; i < fHadronPhys.size(); i++) fHadronPhys[i]->ConstructProcess();

    if (fRestGeant4PhysicsLists->FindPhysicsList("G4RadioactiveDecay")) {
        auto radioactiveDecay = new G4RadioactiveDecay();

        radioactiveDecay->SetHLThreshold(nanosecond);

        // Setting Internal Conversion (ICM) option.
        if (fRestGeant4PhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ICM") == "true")
            radioactiveDecay->SetICM(true);  // Internal Conversion
        else if (fRestGeant4PhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ICM") == "false")
            radioactiveDecay->SetICM(false);  // Internal Conversion
        else if (fRestGeant4PhysicsLists->GetVerboseLevel() >= REST_Essential)
            G4cout << "REST WARNING. restG4. PhysicsList. G4RadioactiveDecay. Option "
                      "ICM not defined."
                   << G4endl;

        // Enabling electron re-arrangment (ARM) option.
        if (fRestGeant4PhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ARM") == "true")
            radioactiveDecay->SetARM(true);  // Internal Conversion
        else if (fRestGeant4PhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ARM") == "false")
            radioactiveDecay->SetARM(false);  // Internal Conversion
        else if (fRestGeant4PhysicsLists->GetVerboseLevel() >= REST_Essential)
            G4cout << "REST WARNING. restG4. PhysicsList. G4RadioactiveDecay. Option "
                      "ARM not defined."
                   << G4endl;
    }

    auto theParticleIterator = GetParticleIterator();

    // to implement UserLimits to StepSize inside the gas
    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4String particleName = particle->GetParticleName();
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
            for (unsigned int n = 0; n < fRestGeant4PhysicsLists->GetIonStepList().size(); n++) {
                if (fRestGeant4PhysicsLists->GetIonStepList()[n] ==
                    G4IonTable::GetIonTable()->GetIonName(Z, A)) {
                    G4ParticleDefinition* particle = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
                    G4String particle_name = G4IonTable::GetIonTable()->GetIonName(Z, A, 0);
                    G4cout << "Found ion: " << particle_name << " Z " << Z << " A " << A << endl;
                    G4ProcessManager* processManager = particle->GetProcessManager();
                    processManager->AddDiscreteProcess(new G4StepLimiter("ionStep"));
                }
            }
        }
}

void PhysicsList::SetCuts() {
    SetDefaultCutValue(1 * mm);
    SetCutsWithDefault();

    spdlog::info("PhysicsList::SetCuts ---> Setting global production cuts from {} to {} keV",
                 fRestGeant4PhysicsLists->GetMinimumEnergyProductionCuts(),
                 fRestGeant4PhysicsLists->GetMaximumEnergyProductionCuts());

    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(
        fRestGeant4PhysicsLists->GetMinimumEnergyProductionCuts() * keV,
        fRestGeant4PhysicsLists->GetMaximumEnergyProductionCuts() * keV);

    SetCutValue(fRestGeant4PhysicsLists->GetCutForGamma() * mm, "gamma");
    SetCutValue(fRestGeant4PhysicsLists->GetCutForElectron() * mm, "e-");
    SetCutValue(fRestGeant4PhysicsLists->GetCutForPositron() * mm, "e+");
    SetCutValue(fRestGeant4PhysicsLists->GetCutForMuon() * mm, "mu+");
    SetCutValue(fRestGeant4PhysicsLists->GetCutForMuon() * mm, "mu-");
    SetCutValue(fRestGeant4PhysicsLists->GetCutForNeutron() * mm, "neutron");
}

void PhysicsList::SetCutValue(G4double cutLengthValue, const G4String& particleName) {
    spdlog::info("PhysicsList::SetCutValue ---> Setting production cut of {:0.2f} mm for particle '{}'",
                 cutLengthValue / mm, particleName);
    G4VUserPhysicsList::SetCutValue(cutLengthValue, particleName);
}
