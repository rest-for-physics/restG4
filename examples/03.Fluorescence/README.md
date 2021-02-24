Date : 24/February/2021

**Detector geometry** : This example implements a simple geometry composed by an empty detector cylinder surrounding the target, which is placed in the center.

**Event generator** : The event generator will launch electrons (electron.rml) or gammas (gamma.rml) in direction to the target. Hits will be registered at the detector.

**Images** directory : contains the spectrum obtained at the detector showing the copper fluorescences.

### Testing the example

```
restG4 gamma.rml (or restG4 electron.rml).
```

```
restManager --c g4Analysis.rml --f Output.root to generate analysis tree entries.
```
