
Int_t Validate(const char* filename) {
    gSystem->Load("libRestFramework.so");
    gSystem->Load("libRestGeant4.so");

    TRestRun run(filename);

    cout << "Entries: " << run.GetEntries() << endl;

    if (run.GetEntries() != 1) {
        cout << "Bad number of entries: " << run.GetEntries() << endl;
        return 1;
    }

    auto metadata = (TRestGeant4Metadata*)run.GetMetadataClass("TRestGeant4Metadata");

    const auto geometryInfo = metadata->GetGeant4GeometryInfo();

    geometryInfo.Print();

    if (geometryInfo.GetAllPhysicalVolumes().size() != 309) {
        cout << "Incorrect number of physical volumes " << geometryInfo.GetAllPhysicalVolumes().size()
             << endl;
        return 2;
    }
    if (geometryInfo.GetAllLogicalVolumes().size() != 22) {
        cout << "Incorrect number of logical volumes " << geometryInfo.GetAllLogicalVolumes().size() << endl;
        return 3;
    }

    if (metadata->GetParticleSource()->GetEnergyDistributionFormulaNPoints() != 2000) {
        cout << "Incorrect number of sampling points on energy distribution: "
             << metadata->GetParticleSource()->GetEnergyDistributionFormulaNPoints() << endl;
        return 4;
    }

    if (metadata->GetParticleSource()->GetAngularDistributionFormulaNPoints() != 300) {
        cout << "Incorrect number of sampling points on angular distribution: "
             << metadata->GetParticleSource()->GetAngularDistributionFormulaNPoints() << endl;
        return 5;
    }

    // event information
    TRestGeant4Event* event = run.GetInputEvent<TRestGeant4Event>();
    run.GetEntry(0);

    if (event->GetNumberOfTracks() != 170) {
        cout << "Incorrect number of tracks: " << event->GetNumberOfTracks() << endl;
        return 6;
    }
    if (event->GetNumberOfHits() != 2933) {
        cout << "Incorrect number of hits: " << event->GetNumberOfHits() << endl;
        return 7;
    }

    constexpr Double_t sensitiveVolumeEnergyRef = 11.4479;
    if (TMath::Abs(event->GetSensitiveVolumeEnergy() - sensitiveVolumeEnergyRef) > 1e-4) {
        cout << "Incorrect sensitive volume energy: " << event->GetSensitiveVolumeEnergy() << endl;
        return 8;
    }

    const auto scintillatorVolumeName =
        "VetoSystem_vetoSystemEast_vetoLayerEast1_assembly-13.veto1_scintillatorVolume-1500.0mm-f1a5df68";
    const auto scintillatorEnergy = event->GetEnergyInVolume(scintillatorVolumeName);
    const auto scintillatorEnergyRef = 2237.1460;
    if (TMath::Abs(scintillatorEnergy - scintillatorEnergyRef) > 1e-4) {
        cout << "Incorrect scintillator volume energy: " << scintillatorEnergy << endl;
        return 9;
    }

    double scintillatorEnergyFromTracks = 0;
    for (const auto& track : event->GetTracks()) {
        const auto hits = track.GetHits();
        for (int i = 0; i < hits.GetNumberOfHits(); i++) {
            const auto volumeName = geometryInfo.GetVolumeFromID(hits.GetVolumeId(i));
            if (volumeName != scintillatorVolumeName) {
                continue;
            }
            const auto energy = hits.GetEnergy(i);
            scintillatorEnergyFromTracks += energy;
        }
    }

    if (TMath::Abs(scintillatorEnergyFromTracks - scintillatorEnergyRef) > 1e-4) {
        cout << "Incorrect scintillator volume energy from tracks: " << scintillatorEnergyFromTracks << endl;
        return 10;
    }

    cout << "All tests passed! [\033[32mOK\033[0m]\n";

    return 0;
}
