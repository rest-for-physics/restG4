//
// Created by lobis on 10/11/2021.
//

#ifndef REST_GLOBALMANAGER_H
#define REST_GLOBALMANAGER_H

#include <TRestGeant4Event.h>
#include <TString.h>

#include <globals.hh>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <utility>
#include <vector>

class TTree;
class TRestGeant4Metadata;
class TRestRun;
class TRestGeant4PhysicsLists;
class TRestGDMLParser;

class GlobalManager {
   public:
    static GlobalManager* Instance();
    ~GlobalManager();

    void InitializeFromConfigFile(const TString&);
    void InitializeRestGeant4Metadata(const TString&);
    void InitializeRestRun(const TString&);
    void InitializeRestGeant4PhysicsLists(const TString&);
    void InitializeTrees();

    inline void SetVolumeLookupTable(std::map<std::string, std::string> table) {
        fVolumeLookupTable = std::move(table);
    }
    inline std::string GetVolumeFromLookupTable(const std::string& key) const {
        if (fVolumeLookupTable.count(key) > 0) {
            return fVolumeLookupTable.at(key);
        }
        return "";
    }

    void WriteGeometry();

    size_t InsertEvent(std::unique_ptr<TRestGeant4Event>&);

    inline TRestGeant4Metadata* GetRestGeant4Metadata() const { return fRestGeant4Metadata; }
    inline TRestRun* GetRestRun() const { return fRestRun; }
    inline TRestGeant4PhysicsLists* GetRestGeant4PhysicsLists() const { return fRestGeant4PhysicsLists; }
    inline TRestGDMLParser* GetRestGDMLParser() const { return fRestGDMLParser; }

    Long64_t GetEntries();

   private:
    GlobalManager();
    static GlobalManager* pinstance_;

    TString fInputConfigFile;

    TRestGeant4Metadata* fRestGeant4Metadata;
    TRestRun* fRestRun;
    TRestGeant4PhysicsLists* fRestGeant4PhysicsLists;
    TRestGDMLParser* const fRestGDMLParser;

    TRestGeant4Event fEvent;
    std::queue<std::unique_ptr<TRestGeant4Event> > fEventContainer;
    std::mutex fEventContainerMutex;

    TTree* fEventTree;
    TTree* fAnalysisTree;

    std::map<std::string, std::string> fVolumeLookupTable;

   public:
    void FillEvents();
    void WriteEvents();

    /*
     * User settings
     */
   private:
    Bool_t fSaveFlag = true;
    Bool_t fSaveAllEventsFlag = false;
    Bool_t fInteractiveFlag = false;

   public:
    inline bool GetSaveFlag() const { return fSaveFlag; }
    inline bool GetSaveAllEventsFlag() const { return fSaveAllEventsFlag; }
    inline bool GetInteractiveFlag() const { return fInteractiveFlag; }

    inline void SetSaveFlag(Bool_t newValue) { fSaveFlag = newValue; }
    inline void SetSaveAllEventsFlag(Bool_t newValue) { fSaveAllEventsFlag = newValue; }
    inline void SetInteractiveFlag(Bool_t newValue) { fInteractiveFlag = newValue; }

    // other
   private:
    TH1D* fPrimaryEnergyDistribution = nullptr;
    TH1D* fPrimaryAngularDistribution = nullptr;
    Double_t fPrimaryEnergyDistributionMin;
    Double_t fPrimaryEnergyDistributionMax;

   public:
    inline TH1D* GetPrimaryEnergyDistribution() const { return fPrimaryEnergyDistribution; }
    inline TH1D* GetPrimaryAngularDistribution() const { return fPrimaryAngularDistribution; }
    inline Double_t GetPrimaryEnergyDistributionMin() const { return fPrimaryEnergyDistributionMin; }
    inline Double_t GetPrimaryEnergyDistributionMax() const { return fPrimaryEnergyDistributionMax; }
};
#endif  // REST_GLOBALMANAGER_H
