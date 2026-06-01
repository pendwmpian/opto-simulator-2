# Wang Fig. 2D-E short-pulse simulation

## Purpose

Use the best Wang Fig. 2A-C voltage-clamp parameter set without refitting `gbar`, then simulate short blue-light pulses as a Fig. 2D-E-style validation of duration dependence.

## Fixed parameters

- Cell: Dura-Bernal PT5B full morphology.
- Biophysics: original Dura-Bernal mechanisms.
- Slice: Wang-like in vitro vertical slice abstraction.
- Temperature: 22 deg C.
- In vitro optical attenuation: `mu_eff = 1.3 mm^-1`.
- ChR2: Foutz/Nikolic 4-state mechanism with Williams Q10 scaling.
- ChR2 density: `gbar = 0.10090350448414476 mS/cm2`.
- Light irradiance: `9.2 mW/mm2`.
- Clamp/readout:
  - SEClamp holding voltage: -70 mV.
  - Rs: 10 MOhm.
  - electrode Rs: 10 MOhm.
  - pipette capacitance: 100 pF.
  - amplifier filter tau: 0.25 ms.
  - amplifier filter order: 4.

This is the same reference set used after Wang Fig. 2A-C calibration. `gbar` is fixed here; it is not fit to the short-pulse results.

## Output

- Metrics: `outputs/wang_fig2de_flash_grid_best_vclamp/fig2de_short_pulse_metrics.csv`
- Summary: `outputs/wang_fig2de_flash_grid_best_vclamp/fig2de_short_pulse_summary.csv`
- Figure: `outputs/wang_fig2de_flash_grid_best_vclamp/fig2de_short_pulse.png`

The plotted current is inverted so inward photocurrent is downward, matching the convention in Wang Fig. 2.

## Result

| duration (ms) | filtered peak (nA) | filtered TTP (ms) | normalized peak | normalized charge |
|---:|---:|---:|---:|---:|
| 1 | 0.149 | 7.25 | 0.256 | 0.138 |
| 2 | 0.343 | 7.73 | 0.589 | 0.323 |
| 3 | 0.460 | 8.13 | 0.791 | 0.441 |
| 4 | 0.522 | 8.48 | 0.898 | 0.383 |
| 5 | 0.553 | 8.75 | 0.951 | 0.412 |
| 6 | 0.569 | 9.03 | 0.978 | 0.431 |
| 7 | 0.577 | 9.30 | 0.992 | 0.445 |
| 8 | 0.581 | 9.55 | 0.998 | 0.457 |
| 10 | 0.582 | 9.72 | 1.000 | 0.478 |
| 20 | 0.582 | 9.72 | 1.000 | 0.565 |
| 50 | 0.582 | 9.72 | 1.000 | 0.754 |
| 100 | 0.582 | 9.72 | 1.000 | 1.000 |

Estimated duration half-response:

- Peak-current half-duration: `2.19 ms`.
- Charge half-duration: `18.49 ms`.

## Interpretation

The model predicts that peak photocurrent saturates quickly with pulse duration. At 9.2 mW/mm2, 5 ms already reaches about 95% of the peak produced by longer pulses, and 8-10 ms reaches the maximum. This is qualitatively compatible with the idea that short flashes can evoke large ChR2 photocurrents once the opening transition is sufficiently driven.

Charge transfer keeps increasing after peak current saturates, because longer light pulses maintain open/desensitized states for longer. Therefore Fig. 2D/E-like traces should be interpreted differently depending on whether the experimental readout is peak current, integrated current, or spike probability.

One caveat is that the simulated peak occurs after light offset for very short pulses. This is expected for this conductance model plus electrode filtering: the photocurrent continues to rise briefly while ChR2 state occupancy and the electrode/amplifier readout catch up. If Wang's short-pulse traces peak substantially earlier or later, that discrepancy would constrain the ChR2/electrode kinetics rather than `gbar`.
