
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
#include <G4PhotonEvaporation.hh>
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
#include <G4UImanager.hh>
#include <G4UnitsTable.hh>
#include <G4UniversalFluctuation.hh>

#include "Particles.h"

using namespace std;

PhysicsList::PhysicsList(SimulationManager* simulationManager, TRestGeant4PhysicsLists* physicsLists)
    : G4VModularPhysicsList(), fSimulationManager(simulationManager) {
    // add new units for radioActive decays
    const G4double G4minute = 60 * second;
    const G4double G4hour = 60 * G4minute;
    const G4double G4day = 24 * G4hour;
    const G4double G4year = 365 * G4day;
    new G4UnitDefinition("minute", "min", "Time", G4minute);
    new G4UnitDefinition("hour", "h", "Time", G4hour);
    new G4UnitDefinition("day", "d", "Time", G4day);
    new G4UnitDefinition("year", "y", "Time", G4year);

    defaultCutValue = 0.1 * mm;

    fRestPhysicsLists = physicsLists;
    G4LossTableManager::Instance();
    // fix lower limit for cut
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(
        fRestPhysicsLists->GetMinimumEnergyProductionCuts() * keV,
        fRestPhysicsLists->GetMaximumEnergyProductionCuts() * keV);

    InitializePhysicsLists();
}

PhysicsList::~PhysicsList() {
    delete fEmPhysicsList;
    delete fDecPhysicsList;
    delete fRadDecPhysicsList;
    for (auto& hadronicPhysicsList : fHadronPhys) {
        delete hadronicPhysicsList;
    }
}

void PhysicsList::InitializePhysicsLists() {
    // Decay physics and all particles
    if (fRestPhysicsLists->FindPhysicsList("G4DecayPhysics") >= 0) {
        fDecPhysicsList = new G4DecayPhysics();
    } else if (fRestPhysicsLists->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        G4cout << "restG4. PhysicsList. G4DecayPhysics is not enabled!!" << G4endl;
    }

    // RadioactiveDecay physicsList
    if (fRestPhysicsLists->FindPhysicsList("G4RadioactiveDecayPhysics") >= 0) {
        fRadDecPhysicsList = new G4RadioactiveDecayPhysics();
    } else if (fRestPhysicsLists->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        G4cout << "restG4. PhysicsList. G4RadioactiveDecayPhysics is not enabled!!" << G4endl;
    }

    // Electromagnetic physicsList
    int emCounter = 0;
    auto emPhysicsListName = "G4EmLivermorePhysics";
    if (fRestPhysicsLists->FindPhysicsList(emPhysicsListName) >= 0) {
        if (fEmPhysicsList == nullptr) {
            fEmPhysicsList = new G4EmLivermorePhysics();
            fEmPhysicsListName = emPhysicsListName;
        }
        emCounter++;
    }

    emPhysicsListName = "G4EmPenelopePhysics";
    if (fRestPhysicsLists->FindPhysicsList(emPhysicsListName) >= 0) {
        if (fEmPhysicsList == nullptr) {
            fEmPhysicsList = new G4EmPenelopePhysics();
            fEmPhysicsListName = emPhysicsListName;
        }
        emCounter++;
    }

    emPhysicsListName = "G4EmStandardPhysics_option3";
    if (fRestPhysicsLists->FindPhysicsList(emPhysicsListName) >= 0) {
        if (fEmPhysicsList == nullptr) {
            fEmPhysicsList = new G4EmStandardPhysics_option3();
            fEmPhysicsListName = emPhysicsListName;
        }
        emCounter++;
    }

    emPhysicsListName = "G4EmStandardPhysics_option4";
    if (fRestPhysicsLists->FindPhysicsList(emPhysicsListName) >= 0) {
        if (fEmPhysicsList == nullptr) {
            fEmPhysicsList = new G4EmStandardPhysics_option4();
            fEmPhysicsListName = emPhysicsListName;
        }
        emCounter++;
    }

    if (fRestPhysicsLists->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Essential &&
        emCounter == 0) {
        RESTWarning << "PhysicsList: No EM physics list has been enabled" << RESTendl;
    }

    if (emCounter > 1) {
        cerr << "PhysicsList: More than 1 EM PhysicsList enabled." << endl;
        exit(1);
    }

    // Hadronic PhysicsList
    if (fRestPhysicsLists->FindPhysicsList("G4HadronPhysicsQGSP_BIC_HP") >= 0) {
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BIC_HP());
    }

    if (fRestPhysicsLists->FindPhysicsList("G4IonBinaryCascadePhysics") >= 0) {
        fHadronPhys.push_back(new G4IonBinaryCascadePhysics());
    }

    if (fRestPhysicsLists->FindPhysicsList("G4HadronElasticPhysicsHP") >= 0) {
        fHadronPhys.push_back(new G4HadronElasticPhysicsHP());
    }

    if (fRestPhysicsLists->FindPhysicsList("G4NeutronTrackingCut") >= 0) {
        fHadronPhys.push_back(new G4NeutronTrackingCut());
    }

    if (fRestPhysicsLists->FindPhysicsList("G4EmExtraPhysics") >= 0) {
        fHadronPhys.push_back(new G4EmExtraPhysics());
    }

    G4cout << "Number of hadronic physics lists added " << fHadronPhys.size() << G4endl;
}

void PhysicsList::ConstructParticle() {
    // pseudo-particles
    G4Geantino::GeantinoDefinition();

    // particles defined in PhysicsLists
    if (fDecPhysicsList) {
        fDecPhysicsList->ConstructParticle();
    }

    if (fEmPhysicsList) {
        fEmPhysicsList->ConstructParticle();
    }

    if (fRadDecPhysicsList) {
        fRadDecPhysicsList->ConstructParticle();
    }

    for (auto& hadronicPhysicsList : fHadronPhys) {
        hadronicPhysicsList->ConstructParticle();
    }
}

void PhysicsList::ConstructProcess() {
    AddTransportation();

    // Electromagnetic physics list
    if (fEmPhysicsList) {
        fEmPhysicsList->ConstructProcess();
        fEmConfig.AddModels();

        G4UImanager* UI = G4UImanager::GetUIpointer();
        UI->ApplyCommand("/process/em/fluo true");
        UI->ApplyCommand("/process/em/auger true");
        UI->ApplyCommand("/process/em/pixe true");

        bool boolEmOptionPixe = StringToBool(
            fRestPhysicsLists->GetPhysicsListOptionValue(fEmPhysicsListName.c_str(), "pixe", "false").Data());
        string stringEmOptionPixe = (boolEmOptionPixe ? "true" : "false");
        G4cout << "Setting EM option '/process/em/pixe' to '" << stringEmOptionPixe << "' for physics list '"
               << fEmPhysicsListName << "'" << endl;
        UI->ApplyCommand(string("/process/em/pixe ") + stringEmOptionPixe);

        bool boolEmOptionFluo = StringToBool(
            fRestPhysicsLists->GetPhysicsListOptionValue(fEmPhysicsListName.c_str(), "fluo", "true").Data());
        string stringEmOptionFluo = (boolEmOptionFluo ? "true" : "false");
        G4cout << "Setting EM option '/process/em/fluo' to '" << stringEmOptionFluo << "' for physics list '"
               << fEmPhysicsListName << "'" << endl;
        UI->ApplyCommand(string("/process/em/fluo ") + stringEmOptionFluo);

        bool boolEmOptionAuger = StringToBool(
            fRestPhysicsLists->GetPhysicsListOptionValue(fEmPhysicsListName.c_str(), "auger", "true").Data());
        string stringEmOptionAuger = (boolEmOptionAuger ? "true" : "false");
        G4cout << "Setting EM option '/process/em/auger' to '" << stringEmOptionAuger
               << "' for physics list '" << fEmPhysicsListName << "'" << endl;
        UI->ApplyCommand(string("/process/em/auger ") + stringEmOptionAuger);
    }

    // Decay physics list
    if (fDecPhysicsList) {
        fDecPhysicsList->ConstructProcess();
    }

    // Radioactive decay
    if (fRadDecPhysicsList) {
        fRadDecPhysicsList->ConstructProcess();
    }

    // Hadronic physics lists
    for (auto& hadronicPhysicsList : fHadronPhys) {
        hadronicPhysicsList->ConstructProcess();
    }

    if (fRestPhysicsLists->FindPhysicsList("G4RadioactiveDecay")) {
        auto radioactiveDecay = new G4RadioactiveDecay();

        const auto decayTimeThreshold = nanosecond;
#ifdef GEANT4_VERSION_LESS_11_0_0
        radioactiveDecay->SetHLThreshold(decayTimeThreshold);
#else
        radioactiveDecay->SetThresholdForVeryLongDecayTime(decayTimeThreshold);
        // ICM option not tested, if it does not work as expected please post an issue on
        // https://github.com/rest-for-physics/restG4/issues
        auto photonEvaporation = new G4PhotonEvaporation();
#endif
        // Setting Internal Conversion (ICM) option.
        if (fRestPhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ICM") == "true") {
#ifdef GEANT4_VERSION_LESS_11_0_0
            radioactiveDecay->SetICM(true);
#else
            photonEvaporation->SetICM(true);
#endif
        } else if (fRestPhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ICM") == "false") {
#ifdef GEANT4_VERSION_LESS_11_0_0
            radioactiveDecay->SetICM(false);
#else
            photonEvaporation->SetICM(false);
#endif
        } else if (fRestPhysicsLists->GetVerboseLevel() >=
                   TRestStringOutput::REST_Verbose_Level::REST_Essential) {
            RESTWarning << "PhysicsList 'G4RadioactiveDecay' option 'ICM' not defined" << RESTendl;
        }

        // Enabling electron re-arrangement (ARM) option.
        if (fRestPhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ARM") == "true") {
            radioactiveDecay->SetARM(true);
        } else if (fRestPhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "ARM") == "false") {
            radioactiveDecay->SetARM(false);
        } else if (fRestPhysicsLists->GetVerboseLevel() >=
                   TRestStringOutput::REST_Verbose_Level::REST_Essential) {
            RESTWarning << "PhysicsList 'G4RadioactiveDecay' option 'ARM' not defined" << RESTendl;
        }

#ifndef GEANT4_VERSION_LESS_11_0_0
        /*
         * If no TritiumDecay is set, do not produce tritium decay unless particle source is H3
         * If H3 is particle source, produce tritium decay unless TritiumDecay option is explicitly disabled
         */
        if (!(fRestPhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "TritiumDecay") ==
              "false") &&
            (fRestPhysicsLists->GetPhysicsListOptionValue("G4RadioactiveDecay", "TritiumDecay", "false") ==
                 "true" ||
             fSimulationManager->GetRestMetadata()->GetParticleSource()->GetParticleName() == "H3")) {
            // Tritium (H3) fix | https://geant4-forum.web.cern.ch/t/triton-decay-with-rdecay01-example/2616
            G4ParticleDefinition* tritium = G4Triton::Definition();
            tritium->SetPDGStable(false);

            G4VProcess* decay = nullptr;
            G4ProcessManager* tritiumProcessManager = tritium->GetProcessManager();
            G4ProcessVector* tritiumProcessVector = tritiumProcessManager->GetAtRestProcessVector();
            for (unsigned int i = 0; i < tritiumProcessVector->size() && decay == nullptr; i++) {
                if ((*tritiumProcessVector)[i]->GetProcessName() == "Decay")
                    decay = (*tritiumProcessVector)[i];
            }
            if (decay) {
                tritiumProcessManager->RemoveProcess(decay);
            }

            decay = new G4RadioactiveDecay();

            tritium->GetProcessManager()->AddProcess(decay, 1000, -1, 1000);
        }
#endif
    }

    auto theParticleIterator = GetParticleIterator();

    // To implement UserLimits to StepSize inside the gas
    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
        G4ParticleDefinition* particle = theParticleIterator->value();
        const auto& particleName = particle->GetParticleName();
        G4ProcessManager* processManager = particle->GetProcessManager();

        if (particleName == "e-") {
            processManager->AddDiscreteProcess(new G4StepLimiter("e-Step"));
        } else if (particleName == "e+") {
            processManager->AddDiscreteProcess(new G4StepLimiter("e+Step"));
        }

        if (particleName == "mu-") {
            processManager->AddDiscreteProcess(new G4StepLimiter("mu-Step"));
        } else if (particleName == "mu+") {
            processManager->AddDiscreteProcess(new G4StepLimiter("mu+Step"));
        }
    }

    // There might be a better way to do this
    for (int Z = 1; Z <= 40; Z++) {
        for (int A = 2 * Z; A <= 3 * Z; A++) {
            for (unsigned int n = 0; n < fRestPhysicsLists->GetIonStepList().size(); n++) {
                if (fRestPhysicsLists->GetIonStepList()[n] == G4IonTable::GetIonTable()->GetIonName(Z, A)) {
                    G4ParticleDefinition* particle = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
                    G4String particle_name = G4IonTable::GetIonTable()->GetIonName(Z, A, 0);
                    cout << "Found ion: " << particle_name << " Z " << Z << " A " << A << endl;
                    G4ProcessManager* processManager = particle->GetProcessManager();
                    processManager->AddDiscreteProcess(new G4StepLimiter("ionStep"));
                }
            }
        }
    }
}

void PhysicsList::SetCuts() {
    SetCutsWithDefault();

    SetCutValue(fRestPhysicsLists->GetCutForGamma() * mm, "gamma");
    SetCutValue(fRestPhysicsLists->GetCutForElectron() * mm, "e-");
    SetCutValue(fRestPhysicsLists->GetCutForPositron() * mm, "e+");
    SetCutValue(fRestPhysicsLists->GetCutForMuon() * mm, "mu+");
    SetCutValue(fRestPhysicsLists->GetCutForMuon() * mm, "mu-");
    SetCutValue(fRestPhysicsLists->GetCutForNeutron() * mm, "neutron");
}
