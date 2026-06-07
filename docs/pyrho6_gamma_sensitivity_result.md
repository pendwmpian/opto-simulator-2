# PyRhO6 gamma sensitivity result

## Purpose

Test whether the low sustained current in PyRhO6 Wang Fig. 2A-C calibration can be explained by the relative conductance of the light-adapted open state `O2`.

This is a sensitivity analysis only. It does not replace the original PyRhO6 fit.

## Setup

- Model: PyRhO6.
- Q10: 1.0 only.
- Gamma scales: `1, 2, 5, 10, 20`.
- Original PyRhO6 `gamma = 0.00369`.
- Effective conductance term: `open = O1 + gamma * O2`.
- For each gamma scale, `gbar` was recalibrated to Wang saturated current target.
- Fig. 2A-C conditions otherwise match the previous PyRhO6 run.

Output:

- `outputs/wang_pyrho6_gamma_sensitivity_q10_1/pyrho6_fig2abc_summary.csv`
- `outputs/wang_pyrho6_gamma_sensitivity_q10_1/pyrho6_fig2abc_intensity.csv`
- `outputs/wang_pyrho6_gamma_sensitivity_q10_1/pyrho6_fig2abc_summary.png`

## Results

| gamma scale | effective gamma | gbar (mS/cm2) | peak 9.2 (nA) | steady 9.2 (nA) | steady/peak | TTP (ms) | intensity K (mW/mm2) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.00369 | 0.1200 | 0.627 | 0.242 | 0.386 | 9.90 | 1.32 |
| 2 | 0.00738 | 0.1200 | 0.629 | 0.249 | 0.395 | 9.90 | 1.31 |
| 5 | 0.01845 | 0.1091 | 0.585 | 0.266 | 0.455 | 10.10 | 1.30 |
| 10 | 0.03690 | 0.1091 | 0.615 | 0.404 | 0.657 | 10.95 | 1.17 |
| 20 | 0.07380 | 0.0677 | 0.613 | 0.589 | 0.960 | 101.30 | 0.65 |

## Interpretation

Increasing `gamma` does directly address the sustained-current mismatch.

- `gamma x10` gives a 1-s steady current around `0.404 nA`, close to the visual Wang Fig. 2A steady-current estimate.
- `gamma x10` keeps time-to-peak near `11 ms`, which is also plausible for Wang Fig. 2A-C.
- `gamma x20` is probably too large: the trace becomes almost plateau-dominated, time-to-peak shifts to about `100 ms`, and the peak/steady distinction mostly disappears.

This means that the sustained-current problem can be explained by the PyRhO6 fitted `O2` conductance being too low for Wang line18 conditions. However, because `gamma` was already fit in the PyRhO source model, this should not be adopted without additional justification.

Practical interpretation:

- Original `gamma` remains the literature baseline.
- `gamma x10` is a useful sensitivity point because it reproduces the Wang-like sustained current.
- If Fig. 2D-F also improves under `gamma x10`, then this suggests Wang line18 behaves as if the light-adapted open state is more conductive than the PyRhO default fit.
- If Fig. 2D-F worsens, the sustained current mismatch should be attributed to another mechanism instead of changing `gamma`.

