# Wang ChR2 Q10 sensitivity result

## Purpose

Test whether the low sustained photocurrent in the Wang Fig. 2 simulation is caused by the temperature correction applied to the ChR2 kinetic rates.

## Fixed setup

- Cell: Dura-Bernal PT5B full morphology.
- Biophysics: original Dura-Bernal mechanisms.
- Temperature: 22 deg C.
- Slice/optics: Wang-like in vitro slice, `mu_eff = 1.3 mm^-1`.
- Voltage clamp/readout: Rs 10 MOhm, electrode Rs 10 MOhm, pipette capacitance 100 pF, amplifier tau 0.25 ms, filter order 4.
- Light: 9.2 mW/mm2 for tests; saturating 100 mW/mm2 for `gbar` calibration.
- `gbar` was re-calibrated separately for each Q10 condition against Wang's reported saturated maximum current.

## Q10 conditions

| condition | meaning |
|---|---|
| none | no 37 -> 22 deg C kinetic scaling |
| generic1.5 | all ChR2 rates scaled by Q10 = 1.5 |
| generic2.0 | all ChR2 rates scaled by Q10 = 2.0 |
| generic2.5 | all ChR2 rates scaled by Q10 = 2.5 |
| williams | parameter-specific Williams-derived scaling |

## Results

| condition | gbar (mS/cm2) | peak 1s (nA) | steady 1s (nA) | steady/peak | TTP (ms) | tau (ms) | peak Kdur (ms) | charge Kdur (ms) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| none | 0.1104 | 0.614 | 0.259 | 0.423 | 8.85 | 101.6 | 1.98 | 36.98 |
| generic1.5 | 0.1009 | 0.607 | 0.238 | 0.392 | 10.70 | 69.4 | 2.38 | 31.66 |
| generic2.0 | 0.0922 | 0.577 | 0.218 | 0.377 | 12.65 | 62.0 | 2.77 | 28.99 |
| generic2.5 | 0.0922 | 0.585 | 0.218 | 0.372 | 14.85 | 68.5 | 3.20 | 27.73 |
| williams | 0.1009 | 0.581 | 0.132 | 0.227 | 9.80 | 27.4 | 2.15 | 18.31 |

Outputs:

- `outputs/wang_q10_sensitivity/q10_sensitivity_summary.csv`
- `outputs/wang_q10_sensitivity/q10_sensitivity_duration.csv`
- `outputs/wang_q10_sensitivity/q10_sensitivity.png`

## Interpretation

The sustained-current problem is strongly Q10-dependent.

The Williams parameter-specific correction gives the smallest sustained current: about 0.13 nA at 1 s, with a steady/peak ratio of 0.23. This is the previous mismatch.

Generic Q10 scaling preserves a much larger sustained component. With generic Q10 = 2.0, time-to-peak becomes close to Wang's reported 12 ms, while the 1-s steady current is about 0.22 nA. With no Q10, steady current rises to 0.26 nA but onset becomes faster.

Therefore the low sustained current is not an unavoidable consequence of the base ChR2 model. It is mainly produced by the parameter-specific Williams correction, especially because it changes the balance between opening, desensitization, and recovery rates rather than simply slowing all rates uniformly.

This does not prove that generic Q10 is correct. It means the Williams scaling should be treated as a major uncertainty for WT ChR2 line18 at 22 deg C. A defensible next step is to keep the literature 37 deg C Foutz/Nikolic parameter set for in vivo simulations, and for Wang in vitro calibration report the 22 deg C `gbar` under both Williams and generic-Q10 assumptions as a sensitivity range.

