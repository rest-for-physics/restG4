[![pipeline status](https://gitlab.cern.ch/rest-for-physics/restg4/badges/master/pipeline.svg)](https://gitlab.cern.ch/rest-for-physics/restg4/-/commits/master)
[![Validation](https://github.com/rest-for-physics/restG4/actions/workflows/validation.yml/badge.svg)](https://github.com/rest-for-physics/restG4/actions/workflows/validation.yml)

## Table of Contents

- Usage
- Structure of the output file

## Usage

Typing `restG4 -h` will show all available options.

### Basic Usage

`restG4 simulation.rml` will launch a simulation using `simulation.rml` as the configuration file. The output file will
be saved to the specified path in the configuration.

### Multithreading

`restG4` makes use of the multithreading capabilities of Geant4 and can be used in multithreading mode to significantly
increase simulation speed.

- `restG4 simulation.rml -t 8` will launch a simulation using 8 worker threads.

By default `restG4` is executed in serial mode (one single thread) and is equivalent to the `-t 0` option.

Simulation results are determined by the random seed and the number of threads, so using different number of threads
with the same seed will produce different results. Using the same seed and same number of threads produces the same
results.

The multithreading implementation of `restG4` is new and although it looks to work properly it may have some bugs. If
you find one of such bugs, please post an issue [here](https://github.com/rest-for-physics/restG4/issues) with
the `multithreading` label.

We encourage our users to use the multithreading feature especially when performing exploratory simulations.

### Overriding values from the RML

The CLI interface allows to set a few different parameters without having to modify the RML.

- `restG4 simulation.rml -o output.root` to specify the output file.
- `restG4 simulation.rml -n 1000` to set the number of processed events
  (`n` in Geant4's `/run/beamOn n` or `nEvents` in the rml configuration).
- `restG4 simulation.rml -g geometry.gdml` to specify the geometry file.

### Ending the simulation early

The simulation run will naturally end when the selected number of events have been processed. They can either be
specified in the RML file or via the `-n` flag on the CLI. However, the user can also specify additional conditions to
end the simulation early.

- `restG4 simulation.rml -N 1000` will signal the simulation to stop once 1000 events have been marked as valid (saved
  on the output file). The final simulation will probably have a slightly higher number of entries than the value
  specified in the case of multithreading as the worker threads may still be working on valid events, or they may not
  have been saved yet. The number of processed events (the equivalent of `nEvents`) will be adjusted to the real number
  of processed events, overriding the selected value of the user, so the resulting simulation would be equivalent to a
  plain simulation using this value as `nEvents`. This holds more accurately the large the `-N` parameter.
- `restG4 simulation.rml -S 1h20m30s` will signal the simulation to stop after 1 hour 20 minutes and 30 seconds have
  elapsed since the start. You can specify any value using this format such as `1h15m`, `30s`, etc.
- Sending the interrupt signal (typically via `CTRL+C`) will signal the simulation to stop, and it will attempt to save
  the results into disk before closing the process.

## Structure of the output file

TODO
