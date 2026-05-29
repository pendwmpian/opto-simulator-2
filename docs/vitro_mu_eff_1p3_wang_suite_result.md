# Wang Suite With In Vitro Optical Attenuation

## Purpose

Re-run the Wang et al. 2007 calibration suite using a more appropriate in vitro optical attenuation value:

```text
mu_eff = 1.3 mm^-1
```

The previous in vitro Wang runs used:

```text
mu_eff = 2.12 mm^-1
```

which was inherited from an in vivo-oriented optical estimate and may include blood/vascular absorption not present in acute slices.

## Model

Fixed settings:

```text
cell: PT5B_full
temperature: 22 degC
ChR2 kinetics: literature 4-state values
recording readout: electrode/amplifier model
SEClamp Rs: 10 Mohm
pipette capacitance: 100 pF
electrode readout tau: 1 ms
amplifier/display filter: 4-stage first-order low-pass, tau = 0.25 ms per stage
gbar calibration: saturating current, 100 mW/mm2, 100 ms, target Imax = 0.642 nA
```

Output:

```text
outputs/wang2007_electrode_mu_eff_1p3/electrode_amplifier_summary.csv
outputs/wang2007_electrode_mu_eff_1p3/electrode_amplifier_summary.png
```

## Result

```text
gbar = 0.13478 mS/cm2

raw SEClamp:
  9.2 mW/mm2 peak = 0.580 nA
  time-to-peak = 5.75 ms
  tau = 22.47 ms
  early peak fraction = 0.312
  intensity K = 0.867 mW/mm2

electrode/amplifier readout:
  9.2 mW/mm2 peak = 0.557 nA
  time-to-peak = 8.05 ms
  tau = 22.43 ms
  early peak fraction = 0.039
  intensity K = 0.762 mW/mm2
  intensity Imax = 0.601 nA
  Hill n = 0.995
  duration K = 1.96 ms
```

## Comparison To Previous `mu_eff = 2.12 mm^-1`

```text
filtered peak:
  2.12 mm^-1: 0.553 nA
  1.30 mm^-1: 0.557 nA

filtered time-to-peak:
  2.12 mm^-1: 8.20 ms
  1.30 mm^-1: 8.05 ms

filtered tau:
  2.12 mm^-1: 22.76 ms
  1.30 mm^-1: 22.43 ms

filtered early peak fraction:
  2.12 mm^-1: 0.036
  1.30 mm^-1: 0.039

filtered intensity K:
  2.12 mm^-1: 0.817 mW/mm2
  1.30 mm^-1: 0.762 mW/mm2

filtered duration K:
  2.12 mm^-1: 2.00 ms
  1.30 mm^-1: 1.96 ms
```

## Interpretation

- Changing `mu_eff` from 2.12 to 1.3 mm^-1 has only a modest effect after `gbar` is refit to the saturating current target.
- The filtered 9.2 mW/mm2 peak current is nearly exactly Wang's 557 pA target.
- Intensity K shifts lower, from 0.817 to 0.762 mW/mm2. This remains close to Wang's 0.84 +/- 0.2 mW/mm2.
- The fast onset component remains suppressed by the electrode/amplifier readout.
- Prolonged kinetics remain the main mismatch:
  - model time-to-peak about 8 ms versus Wang about 12 ms;
  - model tau about 22 ms versus Wang about 48 ms;
  - duration K about 2 ms versus Wang about 3.2 ms.

## Current Recommendation

Use `mu_eff = 1.3 mm^-1` as the default in vitro slice attenuation value going forward, because it is more defensible for acute slices than the previous in vivo-derived value.

Do not use this optical correction to explain the remaining prolonged-kinetics mismatch; it mainly affects current scale and intensity response, not the slow decay kinetics.
