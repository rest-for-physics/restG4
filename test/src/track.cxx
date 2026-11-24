#include <TRestGeant4Event.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4Track.h>
#include <TRestRun.h>
#include <gtest/gtest.h>

#include <G4DynamicParticle.hh>
#include <G4Neutron.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4Transportation.hh>

namespace {

void CheckTrackUpdate(const G4VProcess* process) {
    TRestGeant4Metadata metadata;
    ASSERT_TRUE(metadata.isGlobalTimeReset());
    TRestRun run;
    run.AddMetadata(&metadata);
    TRestGeant4Event event;
    event.InitializeReferences(&run);

    G4Track track(new G4DynamicParticle(G4Neutron::Definition(), G4ThreeVector(0, 0, 1), 2 * MeV), 0,
                  G4ThreeVector());
    track.SetTrackID(1);
    track.SetTrackStatus(fStopAndKill);
    track.AddTrackLength(3 * mm);
    track.SetLocalTime(7 * microsecond);
    G4Step step;
    step.SetTrack(&track);
    track.SetStep(&step);
    step.GetPostStepPoint()->SetProcessDefinedStep(process);

    TRestGeant4Track savedTrack(&track);
    savedTrack.SetEvent(&event);
    savedTrack.UpdateTrack(&track);

    EXPECT_DOUBLE_EQ(savedTrack.GetLength(), 3);
    EXPECT_DOUBLE_EQ(savedTrack.GetTimeLength(), 7);
    EXPECT_EQ(track.GetTrackStatus(), fStopAndKill);
    EXPECT_EQ(savedTrack.GetHits().GetNumberOfHits(), 0);
}

}  // namespace

TEST(restG4, TrackKilledWithoutPostStepProcess) { CheckTrackUpdate(nullptr); }

TEST(restG4, TrackWithNonDecayPostStepProcess) {
    G4Transportation process;
    CheckTrackUpdate(&process);
}
