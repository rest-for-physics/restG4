from REST import ROOT

from pathlib import Path
import numpy as np

files = Path("disk_simple").rglob("*.root")

with open("output.txt", "w") as f:
    f.write("file\twall_time\tsimulation_time\tentries\tprimaries\tsurface\tratio\n")

    for file in files:
        file = str(file)
        print(file)
        run = ROOT.TRestRun(file)
        metadata = run.GetMetadataClass("TRestGeant4Metadata")

        wall_time = metadata.GetSimulationWallTime()
        simulation_time = metadata.GetEquivalentSimulatedTime()
        entries = run.GetEntries()
        primaries = metadata.GetNumberOfEvents()
        surface = metadata.GetGeneratorSurfaceCm2()
        radius = np.sqrt(surface / np.pi)
        reference_radius = 17.320508
        ratio = radius / reference_radius

        print(f"Wall time: {wall_time}")
        print(f"Simulation time: {simulation_time}")
        print(f"Entries: {entries}")
        print(f"Primaries: {primaries}")
        print(f"Surface: {surface}")
        print(f"Ratio: {ratio}")

        f.write(f"{file}\t{wall_time}\t{simulation_time}\t{entries}\t{primaries}\t{surface}\t{ratio}\n")
