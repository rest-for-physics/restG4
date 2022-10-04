Date : 19/July/2019

**Detector geometry** : This example implements a simple geometry composed by a gas volume and a vessel, and a water tank.

**Event generator** : The event generator will launch neutrinoless double beta decays from the gas volume.

**Command** : restG4 NLDBD.rml

**Remark** : The second version of NLDBD.rml is used in order to force loading the remote materials file. Inside the GitHub workflow the remote materials file is downloaded and overwritten (except for the reference pipeline that still needs the old materials file).
