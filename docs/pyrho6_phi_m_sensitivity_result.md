# PyRhO6 phi_m sensitivity result

## Purpose

Test whether the remaining mismatch in the Wang Fig. 2B/C intensity half-max can be explained by PyRhO6 light sensitivity, parameterized by `phi_m`.

This was run after the gamma sensitivity showed that `gamma x5` is a reasonable sustained-current candidate.

## Setup

- Model: PyRhO6.
- Q10: 1.0.
- Gamma scale: 5.0.
- Phi-m scales: `1.0, 0.75, 0.6, 0.5, 0.4`.
- Original `phi_m = 5.02e17 photons/mm2/s`.
- Effective `phi_m = 5.02e17 * phi_m_scale`.
- For each phi-m scale, `gbar` was recalibrated to Wang saturated current target.
- Fig. 2A-C conditions otherwise match the PyRhO6 Wang suite.

Outputs:

- `outputs/wang_pyrho6_phi_m_sensitivity_gamma5_q10_1/pyrho6_fig2abc_summary.csv`
- `outputs/wang_pyrho6_phi_m_sensitivity_gamma5_q10_1/pyrho6_fig2abc_intensity.csv`
- `outputs/wang_pyrho6_phi_m_sensitivity_gamma5_q10_1/pyrho6_fig2abc_summary.png`

## Results

| phi_m scale | effective phi_m | gbar (mS/cm2) | peak 9.2 (nA) | steady 9.2 (nA) | steady/peak | TTP (ms) | intensity K (mW/mm2) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.00 | 5.02e17 | 0.1091 | 0.585 | 0.266 | 0.455 | 10.10 | 1.30 |
| 0.75 | 3.77e17 | 0.1200 | 0.655 | 0.301 | 0.459 | 9.60 | 0.94 |
| 0.60 | 3.01e17 | 0.1200 | 0.665 | 0.308 | 0.463 | 9.30 | 0.74 |
| 0.50 | 2.51e17 | 0.1200 | 0.671 | 0.314 | 0.468 | 9.10 | 0.60 |
| 0.40 | 2.01e17 | 0.1200 | 0.677 | 0.322 | 0.476 | 8.90 | 0.47 |

## Interpretation

Lowering `phi_m` increases light sensitivity as expected.

- `phi_m_scale = 0.75` gives `K = 0.94 mW/mm2` and steady current `~301 pA`.
- `phi_m_scale = 0.60` gives `K = 0.74 mW/mm2` and steady current `~308 pA`.

Therefore the Wang Fig. 2A-C targets are bracketed by `phi_m_scale ~0.6-0.75` under the `Q10=1.0, gamma x5` candidate model.

This is a more interpretable adjustment than changing tissue optics:

- `mu_eff` remains the in vitro slice optical attenuation parameter.
- `phi_m` changes the ChR2 molecule/model light sensitivity.

However, `phi_m` is also a fitted PyRhO parameter, so this should remain a Wang-line18 sensitivity fit unless independently justified by WT ChR2/light-source calibration data.

