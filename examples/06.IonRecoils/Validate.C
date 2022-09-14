#include <TRestGeant4Event.h>

Int_t Validate(const char* filename) {
    cout << "Starting validation for '" << filename << "'" << endl;

    TRestRun run(filename);
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();

    cout << "Testing reading of Geant4 metadata class" << endl;
    TRestGeant4Metadata* geant4Metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");
    if (!geant4Metadata) {
        cout << "Problem reading Geant4 metadata class!" << endl;
        return 6;
    }
    geant4Metadata->PrintMetadata();
    const bool isReferenceGeant4Version = geant4Metadata->GetGeant4Version() == "10.4.3";

    double averageNumberOfHits = 0;
    const double averageNumberOfHitsRef = (!isReferenceGeant4Version) ? 14413.5 : 11223.0;
    const double tolerance = 0.001;

    for (int i = 0; i < run.GetEntries(); i++) {
        run.GetEntry(i);
        averageNumberOfHits += double(event->GetNumberOfHits()) / run.GetEntries();
    }

    cout << "Average number of hits: " << averageNumberOfHits << endl;
    if (TMath::Abs(averageNumberOfHits - averageNumberOfHitsRef) / averageNumberOfHitsRef > tolerance) {
        cout << "The average number of hits does not match the reference value of " << averageNumberOfHitsRef
             << endl;
        return 14;
    }

    Double_t duration = run.GetEndTimestamp() - run.GetStartTimestamp();
    cout << "Duration: " << duration << endl;
    if (duration <= 5 || duration > 200) {
        cout << "The duration of the run in seconds is invalid" << endl;
        return 15;
    }

    cout.precision(12);

    cout << "Time since simulation ended: " << time(NULL) - run.GetEndTimestamp() << endl;
    Double_t delay = time(NULL) - run.GetEndTimestamp();

    if (delay < 0 || delay > 20) {
        cout << "The end timestamp is probably wrong" << endl;
        return 16;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
