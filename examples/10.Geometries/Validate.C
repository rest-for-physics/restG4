
Int_t Validate(const char* filename) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    cout << "Entries: " << run.GetEntries() << endl;

    if (run.GetEntries() != 1000) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 1;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");

    const auto geometryInfo = metadata->GetGeant4GeometryInfo();

    geometryInfo.Print();

    if (geometryInfo.GetAllPhysicalVolumes().size() != 374) {
        cout << "Incorrect number of physical volumes " << geometryInfo.GetAllPhysicalVolumes().size()
             << endl;
        // assembly do not work on older geant4 versions...
        // return 1;
    }
    if (geometryInfo.GetAllLogicalVolumes().size() != 22) {
        cout << "Incorrect number of logical volumes " << geometryInfo.GetAllLogicalVolumes().size() << endl;
        // assembly do not work on older geant4 versions...
        // return 1;
    }

    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    run.GetEntry(0);

    for (const auto& track : event->GetTracks()) {
        const auto hits = track.GetHits();
        for (int i = 0; i < hits.GetNumberOfHits() - 1; i++) {
            const TVector3 position = hits.GetPosition(i);
            const TVector3 nextPosition = hits.GetPosition(i + 1);

            const double time = hits.GetTime(i);
            const double nextTime = hits.GetTime(i + 1);

            const TVector3 momentum = hits.GetMomentumDirection(i);

            const double velocity = (nextPosition - position).Mag() / (nextTime - time);
            constexpr double velocityRef = 299792.458;  // speed of light
            if (TMath::Abs(velocity - velocityRef) / velocityRef > 1e-4) {
                cout << "Incorrect velocity: " << velocity << endl;
                return 2;
            }

            if (TMath::Abs(momentum.Angle(nextPosition - position)) > 1e-4) {
                cout << "Incorrect momentum direction. Angle: " << momentum.Angle(nextPosition - position)
                     << endl;
                return 3;
            }
        }
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";
    return 0;
}
