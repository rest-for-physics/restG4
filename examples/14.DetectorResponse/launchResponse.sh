# Xenon. Following the britannica.com the xenon density is 5.887mg/cm3 at 1atm and 0 degrees celsius
export GDML_TARGET_DENSITY=5.887
# Neon. Following britannica.com the neon density is 0.889mg/cm3 at 1atm and 0 degrees celsius -->
export GDML_QUENCHER_DENSITY=0.8999
export GDML_GAS="XenonNeon"
export GDML_PRESSURE=1.5
export GDML_QUENCHER_PCT=30
export GDML_DRIFT=3
export GDML_DETECTOR_SIZE=6
restG4 DetectorResponse.rml
