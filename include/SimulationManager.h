
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

    TRestRun* fRestRun = nullptr;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists = nullptr;
    TRestGeant4Metadata* fRestGeant4Metadata = nullptr;

    TH1D initialEnergySpectrum;
    TH1D initialAngularDistribution;

    void InitializeOutputManager();
    static OutputManager* GetOutputManager() { return fOutputManager; }

    TRestGeant4Event fEvent;  // Branch on EventTree

    size_t InsertEvent(std::unique_ptr<TRestGeant4Event>& event);

    void WriteEvents();
    void WriteEventsAndCloseFile();

   private:
    static thread_local OutputManager* fOutputManager;
    std::mutex fEventContainerMutex;
    std::queue<std::unique_ptr<TRestGeant4Event> > fEventContainer;
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

   private:
    std::unique_ptr<TRestGeant4Event> fEvent{};
    SimulationManager* fSimulationManager = nullptr;
    std::set<std::string> fActiveVolumes = {};

    friend class StackingAction;
};

#endif  // REST_SIMULATIONMANAGER_H
