
#include "RunAction.h"

#include <PrimaryGeneratorAction.h>
#include <TRestGeant4Metadata.h>
#include <TRestRun.h>

#include <G4PhysicalConstants.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <iomanip>

using namespace std;

extern TRestGeant4Metadata* restG4Metadata;
extern TRestRun* restRun;

RunAction::RunAction(PrimaryGeneratorAction* gen) : G4UserRunAction(), fPrimary(gen) {}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*) {
    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void RunAction::EndOfRunAction(const G4Run* run) {
    G4int nbEvents = run->GetNumberOfEvent();
    if (nbEvents == 0) {
        return;
    }

    G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4String partName = particle->GetParticleName();
    // G4double eprimary = fPrimary->GetParticleGun()->GetParticleEnergy();

    G4cout << "======================== run summary ======================";
    G4cout << "\n" << nbEvents << " Events simulated, " << restRun->GetEntries() << "Events stored\n";
    G4cout << "===========================================================";
    G4cout << G4endl;

    // restG4Metadata->PrintMetadata();

    /*

    G4int prec = 4, wid = prec + 2;
    G4int dfprec = G4cout.precision(prec);

    //particle count
    //
    G4cout << " Nb of generated particles: \n" << G4endl;

    map<G4String,G4int>::iterator it;
    for (it = fParticleCount.begin(); it != fParticleCount.end(); it++) {
        G4String name = it->first;
        G4int count   = it->second;
        G4double eMean = fEmean[name]/count;
        G4double eMin = fEmin[name], eMax = fEmax[name];

        G4cout << "  " << setw(13) << name << ": " << setw(7) << count
            << "  Emean = " << setw(wid) << G4BestUnit(eMean, "Energy")
            << "\t( "  << G4BestUnit(eMin, "Energy")
            << " --> " << G4BestUnit(eMax, "Energy")
            << ")" << G4endl;
    }

    //energy momentum balance
    //

    if (fDecayCount > 0) {
        G4double Ebmean = fEkinTot[0]/fDecayCount;
        G4double Pbmean = fPbalance[0]/fDecayCount;

        G4cout << "\n   Ekin Total (Q): mean = "
            << setw(wid) << G4BestUnit(Ebmean, "Energy")
            << "\t( "  << G4BestUnit(fEkinTot[1], "Energy")
            << " --> " << G4BestUnit(fEkinTot[2], "Energy")
            << ")" << G4endl;

        G4cout << "\n   Momentum balance (excluding gamma desexcitation): mean = "
            << setw(wid) << G4BestUnit(Pbmean, "Energy")
            << "\t( "  << G4BestUnit(fPbalance[1], "Energy")
            << " --> " << G4BestUnit(fPbalance[2], "Energy")
            << ")" << G4endl;
    }

    // remove all contents in fParticleCount
    //
    fParticleCount.clear();
    fEmean.clear();  fEmin.clear(); fEmax.clear();

    // restore default precision
    //
    G4cout.precision(dfprec);
    */
}
