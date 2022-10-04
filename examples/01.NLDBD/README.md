Date : 19/July/2019

**Detector geometry** : This example implements a simple geometry composed by a gas volume and a vessel, and a water tank.

**Event generator** : The event generator will launch neutrinoless double beta decays from the gas volume.

**Command** : restG4 NLDBD.rml

**Remark** : The second version of `NLDBD.rml` (NLDBD2.rml) is used in order to force loading the remote materials file. `NLDBD2.rml` will be used by the `ctests`. Inside the GitHub workflow the remote materials file is downloaded and the local version is overwritten (thus that pipeline simply uses `NLDBD.rml`). The reference pipeline still needs the old materials file, which now is available in the example geometry directory (named `rest.xml`).
