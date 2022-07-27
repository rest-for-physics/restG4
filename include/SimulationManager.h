
#ifndef REST_SIMULATIONMANAGER_H
#define REST_SIMULATIONMANAGER_H

#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PhysicsLists.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>

#include <queue>

class OutputManager;

class SimulationManager {
   public:
    SimulationManager();
    ~SimulationManager();

    void InitializeOutputManager();
    static OutputManager* GetOutputManager() { return fOutputManager; }

    TRestGeant4Event fEvent;  // Branch on EventTree

    size_t InsertEvent(std::unique_ptr<TRestGeant4Event>& event);

    void WriteEvents();
    void WriteEventsAndCloseFile();

   public:
    inline TRestRun* GetRestRun() const { return fRestRun; }
    inline TRestGeant4Metadata* GetRestMetadata() const { return fRestGeant4Metadata; }
    inline TRestGeant4PhysicsLists* GetRestPhysicsLists() const { return fRestGeant4PhysicsLists; }

    inline void SetRestRun(TRestRun* run) { fRestRun = run; }
    inline void SetRestMetadata(TRestGeant4Metadata* metadata) { fRestGeant4Metadata = metadata; }
    inline void SetRestPhysicsLists(TRestGeant4PhysicsLists* physicsLists) {
        fRestGeant4PhysicsLists = physicsLists;
    }

    void EndOfRun();
    inline bool GetAbortFlag() const { return fAbortFlag; }
    void StopSimulation();

    inline double GetElapsedTime() const {
        return 1E-9 * (std::chrono::steady_clock::now().time_since_epoch().count() - fTimeStartUnix);
    }

   private:
    static thread_local OutputManager* fOutputManager;
    std::mutex fSimulationManagerMutex;
    std::queue<std::unique_ptr<TRestGeant4Event> > fEventContainer;

    TRestRun* fRestRun = nullptr;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists = nullptr;
    TRestGeant4Metadata* fRestGeant4Metadata = nullptr;

    int fNumberOfProcessedEvents = 0;
    bool fAbortFlag = false;

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
    OutputManager(const SimulationManager*);
    void UpdateEvent();
    void FinishAndSubmitEvent();

    bool IsEmptyEvent() const;

    bool IsValidEvent() const;
    bool IsValidTrack(const G4Track*) const;  // TODO
    bool IsValidStep(const G4Step*) const;    // TODO

    void RecordTrack(const G4Track*);
    void UpdateTrack(const G4Track*);

    void RecordStep(const G4Step*);

    void AddSensitiveEnergy(Double_t energy, const char* physicalVolumeName);
    void AddEnergyToVolumeForProcess(Double_t energy, const char* volumeName, const char* processName);

    inline bool IsActiveVolume(const char* volumeName) const { return fActiveVolumes.count(volumeName) > 0; }

    inline int GetEventCounter() const { return fProcessedEventsCounter; }

    void BeginOfEventAction();

   private:
    std::unique_ptr<TRestGeant4Event> fEvent{};
    SimulationManager* fSimulationManager = nullptr;
    std::set<std::string> fActiveVolumes = {};

    int fProcessedEventsCounter = 0;

    friend class StackingAction;
};

#endif  // REST_SIMULATIONMANAGER_H
