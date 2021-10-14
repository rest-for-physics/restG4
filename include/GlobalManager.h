//
// Created by lobis on 10/11/2021.
//

#ifndef REST_GLOBALMANAGER_H
#define REST_GLOBALMANAGER_H

#include <TRestGeant4DataEvent.h>
#include <TRestGeant4Event.h>
#include <TString.h>

#include <globals.hh>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
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

    inline void SetVolumeLookupTable(std::map<std::string, std::string> table) { fVolumeLookupTable = table; }
    inline std::string GetVolumeFromLookupTable(std::string key) const {
        if (fVolumeLookupTable.count(key) > 0) {
            return fVolumeLookupTable.at(key);
        }
        return "";
    }

    size_t InsertEvent(std::unique_ptr<TRestGeant4DataEvent>&);

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

    TRestGeant4DataEvent fEvent;
    std::queue<std::unique_ptr<TRestGeant4DataEvent> > fEventContainer;
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
    Bool_t fSaveAllEventsFlag = false;

   public:
    inline bool GetSaveAllEventsFlag() const { return fSaveAllEventsFlag; }
};
#endif  // REST_GLOBALMANAGER_H
