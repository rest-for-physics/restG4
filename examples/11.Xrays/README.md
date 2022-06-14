**Date:** 14/March/2022

**Author:** Javier Galan (javier.galan@unizar.es)

This example uses an x-ray gun with photons in the (0,15)keV range hitting a gas volume filled with Xenon in order to generate an efficiency response matrix.

#### Detector geometry

This example implements a simple geometry based on a cylinder full of a Xenon gas mixture.

#### Testing the example
Execute the following to launch the x-rays

```
restG4 xrays.rml
restManager --c analsis.rml --f RestG4_xrays_00001_ArgonISO.root
```

The script GetQE.C will produce few histograms with the calculated efficiency:

```
restRoot -b -q GetQE.C'("Run_g4Analysis_00001_xrays.root")'
```
