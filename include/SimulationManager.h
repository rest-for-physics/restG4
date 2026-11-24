
#ifndef REST_SIMULATIONMANAGER_H
#define REST_SIMULATIONMANAGER_H

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

#include <atomic>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>

#include "ProcessedEventCounter.h"

class OutputManager;

class SimulationManager {
   public:
    SimulationManager();
    ~SimulationManager();

    void InitializeOutputManager();
    static OutputManager* GetOutputManager() { return fOutputManager; }

    TRestGeant4Event fEvent;  // Branch on EventTree

    void InsertEvent(std::unique_ptr<TRestGeant4Event>& event);

    void WriteEvents();

    std::unique_ptr<std::thread> fPeriodicPrintThread = nullptr;
    std::atomic<bool> fPeriodicPrintThreadEndFlag{false};

   public:
    inline TRestRun* GetRestRun() const { return fRestRun; }
    inline TRestGeant4Metadata* GetRestMetadata() const { return fRestGeant4Metadata; }
    inline TRestGeant4PhysicsLists* GetRestPhysicsLists() const { return fRestGeant4PhysicsLists; }

    inline void SetRestRun(TRestRun* run) { fRestRun = run; }
    inline void SetRestMetadata(TRestGeant4Metadata* metadata) { fRestGeant4Metadata = metadata; }
    inline void SetRestPhysicsLists(TRestGeant4PhysicsLists* physicsLists) {
        fRestGeant4PhysicsLists = physicsLists;
    }

    bool GetPeriodicPrintThreadEndFlag() const {
        return fPeriodicPrintThreadEndFlag.load(std::memory_order_relaxed);
    }

    void BeginOfRunAction();
    void EndOfRunAction();

    inline bool GetAbortFlag() const { return fAbortFlag.load(std::memory_order_relaxed); }
    void StopSimulation();

    inline double GetElapsedTime() const {
        return 1E-9 * (std::chrono::steady_clock::now().time_since_epoch().count() - fTimeStartUnix);
    }

    void SyncStatsFromChild(OutputManager* = fOutputManager);

    int GetNumberOfProcessedEvents() const {
        return fNumberOfProcessedEvents.load(std::memory_order_relaxed);
    }
    int GetNumberOfStoredEvents() const { return fNumberOfStoredEvents.load(std::memory_order_relaxed); }

    std::vector<OutputManager*> GetOutputManagerContainer() const { return fOutputManagerContainer; }

   private:
    static thread_local OutputManager* fOutputManager;
    std::mutex fSimulationManagerMutex;
    std::queue<std::unique_ptr<TRestGeant4Event> > fEventContainer;

    TRestRun* fRestRun = nullptr;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists = nullptr;
    TRestGeant4Metadata* fRestGeant4Metadata = nullptr;

    std::atomic<int> fNumberOfProcessedEvents{0};
    std::atomic<int> fNumberOfStoredEvents{0};

    std::atomic<bool> fAbortFlag{false};

    std::vector<OutputManager*> fOutputManagerContainer = {};
    long fTimeStartUnix = 0;

    /* Primary generation */
   public:
    void InitializeUserDistributions();
    inline const TH1D* GetPrimaryEnergyDistribution() const { return &fPrimaryEnergyDistribution; }
    inline const TH1D* GetPrimaryAngularDistribution() const { return &fPrimaryAngularDistribution; }

   private:
    TH1D fPrimaryEnergyDistribution;
    TH1D fPrimaryAngularDistribution;
};

class OutputManager {
   public:
    struct StepVolumeInfo {
        TString name;
        Int_t id = -1;
        bool active = false;
        bool kill = false;
        bool storageResolved = false;
    };

    OutputManager(const SimulationManager*);
    void UpdateEvent();
    void FinishAndSubmitEvent();

    TRestGeant4Metadata* GetRestMetadata() const { return fSimulationManager->GetRestMetadata(); }
    const StepVolumeInfo& GetStepVolumeInfo(const G4Step* step);

    bool IsEmptyEvent() const;

    bool IsValidEvent() const;
    bool IsValidTrack(const G4Track*) const;  // TODO
    bool IsValidStep(const G4Step*) const;    // TODO

    void RecordTrack(const G4Track*);
    void UpdateTrack(const G4Track*);

    void RecordStep(const G4Step*);

    void AddSensitiveEnergy(Double_t energy);
    void AddEnergyToVolumeForParticleForProcess(Double_t energy, const char* volumeName,
                                                const char* particleName, const char* processName);

    inline int TakeEventCounter() { return fProcessedEventsCounter.Take(); }

    void BeginOfEventAction();

    int GetCurrentEventID() const { return fEvent->GetID(); }

    void SetEventTimeWallPrimaryGeneration(double durationSeconds) {
        fEventTimeWallPrimaryGeneration = durationSeconds;
    }

    void SetEventTimeStart(std::chrono::high_resolution_clock::time_point time) { fEventTimeStart = time; }

   private:
    std::unique_ptr<TRestGeant4Event> fEvent{};
    SimulationManager* fSimulationManager = nullptr;

    // Full paths distinguish daughters shared by multiple placements. Geometry and
    // storage selections remain fixed for this OutputManager's worker lifetime.
    std::unordered_map<std::string, StepVolumeInfo> fStepVolumeCache;

    ProcessedEventCounter fProcessedEventsCounter;

    // start time
    double fEventTimeWallPrimaryGeneration = 0;
    std::chrono::high_resolution_clock::time_point fEventTimeStart;

    void RemoveUnwantedTracks();

    friend class StackingAction;
};

#endif  // REST_SIMULATIONMANAGER_H
