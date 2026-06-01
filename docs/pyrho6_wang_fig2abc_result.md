# PyRhO six-state Wang Fig. 2A-C result

## Scope

This run replaces the Foutz/Nikolic four-state ChR2 mechanism with a PyRhO/Grossman-style six-state ChR2 photocycle and only evaluates Wang Fig. 2A-C-style voltage-clamp behavior.

Fig. 2D-F flash-duration simulations were not run here.

## Assumptions

- ChR2 model: PyRhO/Grossman-style six-state model.
- States: `C1, I1, O1, O2, I2, C2`.
- Conducting occupancy: `open = O1 + gamma * O2`.
- `gbar` is membrane-area normalized maximum conductance in `mS/cm2`; PyRhO's single-channel `g0` is not used.
- PyRhO six-state parameters are used directly as model-native room/fit parameters.
- Light conversion uses photon flux from 470 nm photons:
  - `phi = 1000 * irradiance / E_photon / 1e6`
  - units: photons/mm2/s.
- `E_ChR2 = 0 mV`.
- Q10 is applied as uniform scaling to all kinetic rates, because Williams rate-specific scaling does not map cleanly onto the six-state model.
- Q10 sweep:
  - `Q10 = 1.0`: no scaling.
  - `Q10 = 1.5`: `rate_scale = 0.5443` for 37 -> 22 deg C.
  - `Q10 = 2.0`: `rate_scale = 0.3536` for 37 -> 22 deg C.

## Fixed Wang calibration setup

- Cell: Dura-Bernal PT5B full morphology.
- Biophysics: original Dura-Bernal mechanisms.
- Slice: Wang-like in vitro slice abstraction.
- Optical attenuation: `mu_eff = 1.3 mm^-1`.
- Temperature: 22 deg C.
- Voltage clamp: SEClamp at -70 mV.
- Series resistance: 10 MOhm.
- Electrode/amplifier readout:
  - electrode Rs: 10 MOhm,
  - pipette capacitance: 100 pF,
  - amplifier filter tau: 0.25 ms,
  - filter order: 4.

For each Q10 condition, `gbar` was recalibrated against Wang's saturated maximum current target:

- target `Imax = 0.642 nA`,
- saturating irradiance `100 mW/mm2`,
- duration `100 ms`.

## Outputs

- Script: `scripts/run_wang_pyrho6_fig2abc.py`
- Mechanism: `mod/chr2_pyrho6.mod`
- Summary: `outputs/wang_pyrho6_fig2abc/pyrho6_fig2abc_summary.csv`
- Intensity points: `outputs/wang_pyrho6_fig2abc/pyrho6_fig2abc_intensity.csv`
- Figure: `outputs/wang_pyrho6_fig2abc/pyrho6_fig2abc_summary.png`

## Summary results

| Q10 | rate scale | gbar (mS/cm2) | 9.2 peak (nA) | 9.2 steady (nA) | steady/peak | TTP (ms) | tau (ms) | intensity K (mW/mm2) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.0 | 1.0000 | 0.1158 | 0.6078 | 0.2341 | 0.385 | 9.90 | 74.9 | 1.32 |
| 1.5 | 0.5443 | 0.1078 | 0.5941 | 0.2183 | 0.367 | 13.25 | 59.8 | 1.52 |
| 2.0 | 0.3536 | 0.1028 | 0.5762 | 0.2082 | 0.361 | 17.40 | 72.0 | 1.61 |

## Interpretation

Compared with the previous Williams-scaled four-state model, PyRhO6 increases the 1-s sustained current from about 0.13 nA to about 0.21-0.23 nA under these assumptions. This supports the idea that the Foutz/Nikolic + Williams combination was over-desensitizing the model.

However, PyRhO6 still does not recover the visually estimated Wang Fig. 2A steady current around 0.4 nA. It also predicts intensity half-max values around `1.3-1.6 mW/mm2`, higher than Wang's reported/estimated `~0.84 mW/mm2`.

The Q10 tradeoff is clear:

- `Q10 = 1.0` gives the largest sustained current and fastest onset.
- `Q10 = 1.5` gives TTP near the Wang 12 ms range, but sustained current remains low.
- `Q10 = 2.0` slows onset too much and slightly lowers sustained current.

Current best interpretation:

1. PyRhO6 is a better baseline than the previous four-state Williams-scaled model for sustained current.
2. The remaining mismatch likely requires either:
   - different six-state parameter set for WT ChR2 / line18,
   - different conductance ratio `gamma`,
   - different high-light adaptation balance,
   - or a joint constraint using Fig. 2D-F after Fig. 2A-C is accepted.
3. For now, `Q10 = 1.0` and `Q10 = 1.5` bracket the most plausible Wang 22 deg C behavior; `Q10 = 2.0` is useful as a slow-kinetics sensitivity bound.

