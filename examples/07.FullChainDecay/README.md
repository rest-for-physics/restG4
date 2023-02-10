**Date:** 12/March/2021

**Author:** Javier Galan (javier.galan@unizar.es)

This example tests a radiative chain of a full decay chain of U238. Its daughters, and the correct separation into subevents.

#### Detector geometry

This example implements a simple geometry based on a box made of lead.

#### Event generator
The event generator will produce isotope decays from the center of the box. The parent isotope may be changed through the REST_ISOTOPE variable that takes U238 by default.

#### Testing the example
Execute the following to launch a full chain decay till the chain reaches the stable element.

```
restG4 fullChain.rml
```

Execute the following to launch a single decay of the isotope defined by REST_ISOTOPE

```
restG4 singleDecay.rml
```

See also `Validation.C` for validation routines executed in the pipeline.
