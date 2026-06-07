# Wang 4-state 22C no-truncate rerun result

## Purpose

Rerun the Wang Fig. 2A-C voltage-clamp calibration after removing the erroneous `truncate_radius_um=500` morphology cutoff.

Two illumination conditions were tested:

1. `full`: all retained non-axon morphology receives the large-field in vitro illumination.
2. `apical_0p4mm2`: a circular `0.4 mm2` field is placed on the slice face over the apical dendrite area centroid.

Fig. 2D-F was not run.

## Fixed assumptions

- Cell: Dura-Bernal PT5B full morphology.
- Axon excluded from ChR2 rows, as before.
- Slice: 300 um Wang-like in vitro slice geometry.
- `truncate_radius_um`: **none**.
- Temperature: 22 deg C.
- ChR2: previous 4-state mechanism.
- Q10: Williams parameter-specific 37 -> 22 deg C scaling, same as the previous 4-state Wang suite.
- Optical attenuation: `mu_eff = 1.3 mm^-1`.
- Clamp/readout: SEClamp at -70 mV, Rs 10 MOhm, electrode Rs 10 MOhm, pipette capacitance 100 pF, amplifier tau 0.25 ms, filter order 4.
- Calibration target: Wang saturated maximum current `Imax = 0.642 nA`, using 100 mW/mm2 for 100 ms.

## Output

- Script: `scripts/rerun_wang_4state_no_truncate.py`
- Summary: `outputs/wang_4state_22c_no_truncate_rerun/summary.csv`
- Intensity points: `outputs/wang_4state_22c_no_truncate_rerun/intensity.csv`
- Figure: `outputs/wang_4state_22c_no_truncate_rerun/summary.png`

## Area sanity check

| mode | retained area (um2) | illuminated area (um2) | illuminated fraction | apic area (um2) | illuminated apic fraction |
|---|---:|---:|---:|---:|---:|
| full | 26458.9 | 26458.9 | 1.000 | 15432.8 | 1.000 |
| apical_0p4mm2 | 26458.9 | 15621.8 | 0.590 | 15432.8 | 0.970 |

This confirms that the morphology is no longer being truncated by an arbitrary radius. The `apical_0p4mm2` field covers almost all apical dendrite area while excluding part of basal/somatic morphology.

## Results

| mode | gbar (mS/cm2) | peak 9.2 (nA) | steady 9.2 (nA) | steady/peak | TTP (ms) | tau (ms) | K (mW/mm2) |
|---|---:|---:|---:|---:|---:|---:|---:|
| full | 0.1044 | 0.599 | 0.136 | 0.227 | 9.80 | 27.4 | 0.860 |
| apical_0p4mm2 | 0.4401 | 0.589 | 0.136 | 0.230 | 10.15 | 75.6 | 0.743 |

## Interpretation

Removing the 500 um radius cutoff fixes the morphology issue and substantially changes the calibration context:

- retained area increases to about `26459 um2`, compared with about `19910 um2` in the erroneous truncated setup;
- full-field `gbar` is `0.104 mS/cm2`, close to the previous apparent value because the calibration target is peak-dominated;
- the 4-state Williams-scaled model still produces too little 1-s sustained current, around `136 pA`;
- `full` illumination gives an intensity half-max close to Wang's `~0.84 mW/mm2`;
- apical 0.4 mm2 illumination shifts the fitted `gbar` upward because less total membrane area is illuminated, but does not fix the sustained-current mismatch.

Therefore the earlier morphology cutoff was a serious modeling error, but the 4-state sustained-current problem remains even after restoring the morphology.

