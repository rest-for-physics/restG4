
#include "SimulationManager.h"

#include <G4EventManager.hh>
#include <G4Nucleus.hh>
#include <G4Threading.hh>
#include <Randomize.hh>
#include <csignal>

#include "SteppingAction.h"

using namespace std;

thread_local OutputManager* SimulationManager::fOutputManager = nullptr;

SimulationManager::SimulationManager() {
    // Only master thread should create the SimulationManager
    if (!G4Threading::IsMasterThread()) {
        cout << "Only master thread should create the SimulationManager!" << endl;
        exit(1);
    }
}

void SimulationManager::InitializeOutputManager() {
    lock_guard<std::mutex> guard(fSimulationManagerMutex);
    fOutputManager = new OutputManager(this);
    fOutputManagerContainer.push_back(fOutputManager);
}

void PeriodicPrint(SimulationManager* simulationManager) {
    const auto restG4Metadata = simulationManager->GetRestMetadata();

    while (!simulationManager->GetPeriodicPrintThreadEndFlag()) {
        std::this_thread::sleep_for(std::chrono::seconds(2));

        for (auto& outputManager : simulationManager->GetOutputManagerContainer()) {
            simulationManager->SyncStatsFromChild(outputManager);
        }

        string outputPercentageType = "events";
        double completionPercentage = double(simulationManager->GetNumberOfProcessedEvents()) /
                                      double(restG4Metadata->GetNumberOfEvents()) * 100;
        if (restG4Metadata->GetNumberOfRequestedEntries() > 0) {
            auto completionPercentageEntries = double(simulationManager->GetNumberOfStoredEvents()) /
                                               double(restG4Metadata->GetNumberOfRequestedEntries()) * 100;
            if (completionPercentageEntries > completionPercentage) {
                completionPercentage = completionPercentageEntries;
                outputPercentageType = "entries";
            }
        }
        if (restG4Metadata->GetSimulationMaxTimeSeconds() > 0) {
            auto completionPercentageTime = double(simulationManager->GetElapsedTime()) /
                                            double(restG4Metadata->GetSimulationMaxTimeSeconds()) * 100;
            if (completionPercentageTime > completionPercentage) {
                completionPercentage = completionPercentageTime;
                outputPercentageType = "time";
            }
        }

        string progress;
        if (outputPercentageType == "events") {
            progress = " | " + to_string(simulationManager->GetNumberOfProcessedEvents()) +
                       " events processed / " + to_string(restG4Metadata->GetNumberOfEvents()) +
                       " requested (" +
                       TString::Format("%.2e", simulationManager->GetNumberOfProcessedEvents() /
                                                   simulationManager->GetElapsedTime())
                           .Data() +
                       "/s)";
        } else if (outputPercentageType == "entries") {
            progress = " | " + to_string(simulationManager->GetNumberOfStoredEvents()) + " entries / " +
                       to_string(restG4Metadata->GetNumberOfRequestedEntries()) + " requested";
        } else if (outputPercentageType == "time") {
            progress = " | " + ToTimeStringLong(simulationManager->GetElapsedTime()) + " elapsed / " +
                       ToTimeStringLong(simulationManager->GetRestMetadata()->GetSimulationMaxTimeSeconds());
        }

        string timeInfo;
        if (outputPercentageType != "time") {
            timeInfo = " | " + ToTimeStringLong(simulationManager->GetElapsedTime()) + " elapsed";
        }

        G4cout << TString::Format("%5.2f", completionPercentage).Data() << "% | "
               << simulationManager->GetNumberOfStoredEvents() << " events stored / "
               << simulationManager->GetNumberOfProcessedEvents() << " processed ("
               << TString::Format("%.2e", simulationManager->GetNumberOfStoredEvents() /
                                              simulationManager->GetElapsedTime())
                      .Data()
               << "/s)" << progress << timeInfo << G4endl;
    }
}

int interruptSignalHandler(const int, void* ptr) {
    // See https://stackoverflow.com/a/43400143/11776908
    cout << "Stopping Run! Program was manually stopped by user (CTRL+C)!" << endl;
    const auto manager = (SimulationManager*)(ptr);
    manager->StopSimulation();
    return 0;
}

void SimulationManager::BeginOfRunAction() {
    if (G4Threading::IsMultithreadedApplication() && G4Threading::G4GetThreadId() != -1) {
        return;  // Only call this once from the main thread
    }
    signal(SIGINT, (void (*)(int))interruptSignalHandler);  // Add custom signal handler before simulation

    fTimeStartUnix = chrono::steady_clock::now().time_since_epoch().count();
#ifndef GEANT4_WITHOUT_G4RunManagerFactory
    // gives segfault in old Geant4 versions such as 10.4.3, didn't look into it
    if (GetRestMetadata()->PrintProgress() ||
        GetRestMetadata()->GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Essential) {
        fPeriodicPrintThread = make_unique<thread>(&PeriodicPrint, this);
    }
#endif
}

void SimulationManager::EndOfRunAction() {
    if (G4Threading::IsMultithreadedApplication() && G4Threading::G4GetThreadId() != -1) {
        return;  // Only call this once from the main thread
    }

    fPeriodicPrintThreadEndFlag = true;

    WriteEvents();

    if (fPeriodicPrintThread != nullptr) {
        if (fPeriodicPrintThread->joinable()) {
            fPeriodicPrintThread->join();  // need to join thread, it may block for up to 1 thread period
        }
    }

    for (auto& outputManager : fOutputManagerContainer) {
        fNumberOfProcessedEvents += outputManager->GetEventCounter();
        delete outputManager;
    }
    GetRestMetadata()->SetNumberOfEvents(fNumberOfProcessedEvents);

    fOutputManagerContainer.clear();
}

SimulationManager::~SimulationManager() {
    delete fRestRun;
    delete fRestGeant4Metadata;
    delete fRestGeant4PhysicsLists;

    for (auto& outputManager : fOutputManagerContainer) {
        delete outputManager;
    }
}

void SimulationManager::InsertEvent(std::unique_ptr<TRestGeant4Event>& event) {
    lock_guard<mutex> guard(fSimulationManagerMutex);
    fEventContainer.push(std::move(event));
}

void SimulationManager::WriteEvents() {
    lock_guard<mutex> guard(fSimulationManagerMutex);

    if (fEventContainer.empty()) {
        return;
    }

    while (!fEventContainer.empty()) {
        fEvent = *fEventContainer.front();

        if (!fRestGeant4Metadata->GetStoreTracks()) {
            fEvent.ClearTracks();
        }

        const auto eventTree = fRestRun->GetEventTree();
        if (eventTree != nullptr) {
            eventTree->Fill();
        }

        const auto analysisTree = fRestRun->GetAnalysisTree();
        if (analysisTree != nullptr) {
            analysisTree->SetEventInfo(&fEvent);
            analysisTree->Fill();
        }

        fEventContainer.pop();
    }

    const auto nRequestedEntries = GetRestMetadata()->GetNumberOfRequestedEntries();
    if (nRequestedEntries > 0 && !fAbortFlag && fRestRun->GetEventTree()->GetEntries() >= nRequestedEntries) {
        G4cout << "Stopping Run! We have reached the number of requested entries (" << nRequestedEntries
               << ")" << endl;
        StopSimulation();
    }
}

void SimulationManager::InitializeUserDistributions() {
    auto random = []() { return (double)G4UniformRand(); };

    for (int i = 0; i < fRestGeant4Metadata->GetNumberOfSources(); i++) {
        fRestGeant4Metadata->GetParticleSource(i)->SetRandomMethod(random);
    }

    TRestGeant4ParticleSource* source = fRestGeant4Metadata->GetParticleSource(0);

    if (TRestGeant4PrimaryGeneratorTypes::StringToEnergyDistributionTypes(
            source->GetEnergyDistributionType().Data()) ==
        TRestGeant4PrimaryGeneratorTypes::EnergyDistributionTypes::TH1D) {
        TFile file(source->GetEnergyDistributionFilename());
        auto distribution = (TH1D*)file.Get(source->GetEnergyDistributionNameInFile());

        if (!distribution) {
            RESTError << "Error when trying to find energy spectrum" << RESTendl;
            RESTError << "File: " << source->GetEnergyDistributionFilename() << RESTendl;
            RESTError << "Spectrum name: " << source->GetEnergyDistributionNameInFile() << RESTendl;
            exit(1);
        }

        fPrimaryEnergyDistribution = *distribution;
    }

    if (TRestGeant4PrimaryGeneratorTypes::StringToAngularDistributionTypes(
            source->GetAngularDistributionType().Data()) ==
        TRestGeant4PrimaryGeneratorTypes::AngularDistributionTypes::TH1D) {
        TFile file(source->GetAngularDistributionFilename());
        auto distribution = (TH1D*)file.Get(source->GetAngularDistributionNameInFile());

        if (!distribution) {
            RESTError << "Error when trying to find angular spectrum" << RESTendl;
            RESTError << "File: " << source->GetAngularDistributionFilename() << RESTendl;
            RESTError << "Spectrum name: " << source->GetAngularDistributionNameInFile() << RESTendl;
            exit(1);
        }

        fPrimaryAngularDistribution = *distribution;
    }
}

void SimulationManager::StopSimulation() {
    // Still needs to be propagated to other threads, this is done in the BeginOfEventAction
    G4RunManager::GetRunManager()->AbortRun(true);
    fAbortFlag = true;
}

void SimulationManager::SyncStatsFromChild(OutputManager* outputManager) {
    if (outputManager == nullptr) {
        return;
    }
    lock_guard<mutex> guard(fSimulationManagerMutex);
    fNumberOfProcessedEvents += outputManager->GetEventCounter();
    outputManager->ResetEventCounter();
    fNumberOfStoredEvents = fRestRun->GetEntries();
}

// OutputManager
OutputManager::OutputManager(const SimulationManager* simulationManager)
    : fSimulationManager(const_cast<SimulationManager*>(simulationManager)) {
    // this class should only exist on the threads performing the simulation
    if (G4Threading::IsMasterThread() && G4Threading::IsMultithreadedApplication()) {
        G4cout << "Error in 'OutputManager', this instance should never exist" << endl;
        exit(1);
    }
}

void OutputManager::BeginOfEventAction() {
    // This should only be executed once at BeginOfEventAction
    UpdateEvent();
    fProcessedEventsCounter++;

    if (fSimulationManager->GetAbortFlag()) {
        G4RunManager::GetRunManager()->AbortRun(true);
    }

    if (fSimulationManager->GetRestMetadata()->GetSimulationMaxTimeSeconds() != 0 &&
        !fSimulationManager->GetAbortFlag() &&
        fSimulationManager->GetElapsedTime() >
            fSimulationManager->GetRestMetadata()->GetSimulationMaxTimeSeconds()) {
        G4cout << "Stopping Run! We have reached the time limit of "
               << ToTimeStringLong(fSimulationManager->GetRestMetadata()->GetSimulationMaxTimeSeconds())
               << endl;
        fSimulationManager->StopSimulation();
    }
}

void OutputManager::UpdateEvent() {
    auto event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    fEvent = make_unique<TRestGeant4Event>(event);
    fEvent->InitializeReferences(fSimulationManager->GetRestRun());

    if ((event->GetPrimaryVertex()->GetNumberOfParticle() !=
         fEvent->GetGeant4Metadata()->GetNumberOfSources()) &&
        (int(fEvent->GetNumberOfPrimaries()) != fEvent->GetGeant4Metadata()->GetNumberOfSources())) {
        cout << "TRestGeant4Event: Number of particles on primary vertex does not match number of sources "
                "for event ID: "
             << fEvent->GetID() << endl;
        exit(1);
    }

    for (const TString& volumeName : fEvent->GetGeant4Metadata()->GetActiveVolumes()) {
        fEvent->AddActiveVolume(volumeName.Data());
    }
}

bool OutputManager::IsEmptyEvent() const { return !fEvent || fEvent->fTracks.empty(); }

bool OutputManager::IsValidEvent() const {
    if (IsEmptyEvent()) {
        return false;
    }
    if (fSimulationManager->GetRestMetadata()->GetSaveAllEvents()) {
        return true;
    }
    const auto energy = fEvent->GetSensitiveVolumeEnergy();
    if (energy <= 0) {
        return false;
    }
    // TODO: additional checks
    if (energy < fSimulationManager->GetRestMetadata()->GetMinimumEnergyStored() ||
        energy > fSimulationManager->GetRestMetadata()->GetMaximumEnergyStored()) {
        return false;
    }
    return true;
}

void OutputManager::FinishAndSubmitEvent() {
    if (IsValidEvent()) {
        for (int i = 0; i < fEvent->GetNumberOfActiveVolumes(); i++) {
            fEvent->SetEnergyDepositedInVolume(i, fEvent->GetEnergyInVolume(fEvent->fVolumeStoredNames[i]));
        }
        if (fSimulationManager->GetRestMetadata()->GetRemoveUnwantedTracks()) {
            RemoveUnwantedTracks();
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - fEventTimeStart;

        fEvent->fEventTimeWall = elapsed.count() / 1000;  // seconds
        fEvent->fEventTimeWallPrimaryGeneration = fEventTimeWallPrimaryGeneration;

        fSimulationManager->InsertEvent(fEvent);
        fSimulationManager->WriteEvents();
    }
    UpdateEvent();
}

void OutputManager::RecordTrack(const G4Track* track) {
    if (!IsValidTrack(track)) {
        return;
    }
    fEvent->InsertTrack(track);

    if (fEvent->fSubEventID > 0) {
        //       const auto& lastTrack = ;
        assert(fEvent->fTracks.back().GetTrackID() == track->GetTrackID());
        // TODO
        /*
        bool isSubEventPrimary = fEvent->IsTrackSubEventPrimary(lastTrack.fTrackID);
        if (isSubEventPrimary) {
            spdlog::debug(
                "OutputManager::RecordTrack - Setting track ID {} as SubEventPrimaryTrack of EventID {} "
                "(SubEventID {}). Track info: {} - Created by "
                "{} - ParentID: {}",
                lastTrack.fTrackID, fEvent->fEventID, fEvent->fSubEventID, lastTrack.fParticleName,
                lastTrack.fCreatorProcess, lastTrack.fParentID);
        }
         */
    }
}

void OutputManager::UpdateTrack(const G4Track* track) {
    if (!IsValidTrack(track)) {
        return;
    }
    fEvent->UpdateTrack(track);
}

void OutputManager::RecordStep(const G4Step* step) { fEvent->InsertStep(step); }

void OutputManager::AddSensitiveEnergy(Double_t energy) {
    fEvent->AddEnergyToSensitiveVolume(energy);
    /*
        const TString physicalVolumeNameNew = fSimulationManager->GetRestMetadata()->GetGeant4GeometryInfo()
                                                  .GetAlternativeNameFromGeant4PhysicalName(physicalVolumeName);
                                                  */
}

void OutputManager::AddEnergyToVolumeForParticleForProcess(Double_t energy, const char* volumeName,
                                                           const char* particleName,
                                                           const char* processName) {
    fEvent->AddEnergyInVolumeForParticleForProcess(energy, volumeName, particleName, processName);
}
