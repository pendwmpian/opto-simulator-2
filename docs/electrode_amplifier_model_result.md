# Electrode / Amplifier Readout Model

## Purpose

The Wang-suite ChR2 model produced a small fast current immediately after light onset. This may be a real ChR2 `C1 -> O1` fast component that is not visible in Wang et al. 2007 because of recording-chain filtering.

This test adds a simple explicit recording readout model before changing ChR2 kinetics.

## Model

This is not yet a full pipette-electrode feedback circuit. It keeps the somatic `SEClamp` as the voltage-clamp actuator, but makes two recording-chain changes:

```text
SEClamp series/access resistance:
  Rs = 10 Mohm

pipette/electrode RC readout:
  R_electrode = 10 Mohm
  C_pipette = 100 pF
  tau = R * C = 1.0 ms

amplifier/display filter:
  4 cascaded first-order low-pass stages
  tau = 0.25 ms per stage
```

The intent is to approximate pipette capacitance, access filtering, amplifier low-pass behavior, and figure-level smoothing without changing molecular ChR2 parameters.

## Protocol

```text
cell:
  PT5B_full

temperature:
  22 degC

ChR2:
  Foutz/Nikolic baseline rates
  Williams Q10 mapping
  photon-flux activation
  uniform membrane localization

gbar:
  fit to Wang saturating current target
  100 mW/mm2, 100 ms, target 0.642 nA

validation:
  9.2 mW/mm2, 1 s
  Wang 7-point intensity response
  duration-response series
```

Output:

```text
outputs/wang2007_electrode_amplifier_model/electrode_amplifier_summary.csv
outputs/wang2007_electrode_amplifier_model/electrode_amplifier_summary.png
```

## Result

```text
fitted gbar:
  0.13478 mS/cm2

raw SEClamp, before electrode/amplifier readout:
  9.2 mW/mm2 peak = 0.574 nA
  time-to-peak = 5.9 ms
  inactivation tau = 22.8 ms
  early peak fraction = 0.295
  intensity K = 0.923 mW/mm2
  intensity Imax = 0.633 nA
  Hill n = 0.963

filtered electrode/amplifier readout:
  9.2 mW/mm2 peak = 0.553 nA
  time-to-peak = 8.2 ms
  inactivation tau = 22.8 ms
  early peak fraction = 0.036
  intensity K = 0.817 mW/mm2
  intensity Imax = 0.599 nA
  Hill n = 1.00
  duration K = 2.00 ms
```

## Interpretation

- The recording-chain model strongly suppresses the immediate fast component:
  - early peak fraction drops from about 0.30 to about 0.04.
- The filtered 9.2 mW/mm2 peak current is 0.553 nA, very close to Wang's 0.557 nA mean.
- The filtered intensity K is 0.817 mW/mm2, close to Wang's 0.84 +/- 0.2 mW/mm2.
- The apparent time-to-peak moves from about 5.9 ms to 8.2 ms. This is closer to Wang but still below the reported prolonged-stimulation value of about 12 ms.
- The inactivation tau remains too fast at about 22.8 ms, so recording-chain filtering alone does not explain the full prolonged-trace kinetics.
- The duration K remains low at about 2.0 ms versus Wang's about 3.2 ms.

## Current Conclusion

Recording-chain filtering is a strong candidate explanation for why the fast onset component is not visible in Wang's figure. It also brings the current scale and intensity K into good agreement without changing ChR2 kinetics.

The remaining mismatch is mainly the slower prolonged-stimulation kinetics:

```text
time-to-peak:
  model filtered = 8.2 ms
  Wang reported = about 12 ms

inactivation tau:
  model filtered = 22.8 ms
  Wang reported = about 48 ms

duration K:
  model filtered = 2.0 ms
  Wang reported = about 3.2 ms
```

The next refinement, if needed, should target prolonged desensitization/recovery kinetics or a more complete electrode-feedback circuit, not the initial fast component.

## Original-Biophysics Reference Result

The current reference condition uses the original Dura-Bernal PT5B active mechanisms, applies section-level ion parameters from `cellParams.pkl`, subtracts a matched no-light baseline trace, and keeps the same Wang-style recording readout:

```text
SEClamp series/access resistance:
  Rs = 10 Mohm

pipette/electrode RC readout:
  R_electrode = 10 Mohm
  C_pipette = 100 pF
  tau = 1.0 ms

amplifier/display filter:
  4 cascaded first-order low-pass stages
  tau = 0.25 ms per stage

optical condition:
  in vitro mu_eff = 1.3 mm^-1
```

Output:

```text
outputs/wang2007_electrode_amplifier_original_biophys_Rs10_tau1/electrode_amplifier_summary.csv
outputs/wang2007_electrode_amplifier_original_biophys_Rs10_tau1/electrode_amplifier_summary.png
```

Result:

```text
fitted gbar:
  0.10090 mS/cm2

raw SEClamp, before electrode/amplifier readout:
  9.2 mW/mm2 peak = 0.590 nA
  time-to-peak = 7.65 ms
  inactivation tau = 27.19 ms
  early peak fraction = 0.278
  intensity K = 0.906 mW/mm2

filtered electrode/amplifier readout:
  9.2 mW/mm2 peak = 0.581 nA
  time-to-peak = 9.80 ms
  inactivation tau = 27.39 ms
  early peak fraction = 0.041
  intensity K = 0.864 mW/mm2
  duration K = 2.15 ms

voltage-clamp quality:
  soma voltage range = -71.29 to -65.60 mV
```

This is now the preferred Wang Fig. 2A-C reference value for the in vitro PT5B calibration. The voltage clamp is imperfect, but that is expected for this access resistance and large distributed photocurrent. The key point is that using a realistic recording readout moves the apparent intensity K close to Wang's reported 0.84 +/- 0.2 mW/mm2 while preserving original active conductances.
