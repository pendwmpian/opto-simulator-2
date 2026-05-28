# In Vitro Simulation TODO

## Current Scope

The in vitro branch is used to reproduce Wang et al. 2007 line 18 acute-slice photocurrent data and estimate an effective ChR2 `gbar` before returning to in vivo patterned stimulation.

## TODO

- Keep Wang et al. 2007 as the calibration target: acute cortical parasagittal slices, 250-350 um thickness, 465-495 nm light, 9.2 mW/mm2, 1 s pulse, -70 mV voltage clamp.
- Treat slice geometry as a vertical cortical slice containing the cortical depth axis and one horizontal axis. In current morphology coordinates, use x-y as the slice face and z as the 300 um thickness axis.
- Keep the current large-field illumination approximation for now. Wang used large illuminated areas around 0.4 mm2; this likely covers most of the retained single-cell morphology, so exact field-stop geometry can be deferred.
- Add explicit finite 0.4 mm2 circular or square illumination later if calibration becomes sensitive to distal dendrite inclusion.
- Sweep soma depth from the light-entry surface. Acute-slice patching usually avoids cells immediately at the cut surface; 50-100 um below the surface is a common practical range for healthier visually targeted cells.
- Use the soma-depth sweep to choose a representative in vitro calibration position for future Wang-style `gbar` fitting.
- After geometry and current scale are stable, tune ChR2 kinetics to match Wang time-to-peak and inactivation time constant.
- After kinetics are tuned, refit `gbar` and promote the final value to the default `line18_wang2007` expression scale.
- Implement a Wang protocol suite rather than fitting a single trace:
  - Fig. 2A/text scale: 1 s, 9.2 mW/mm2, peak photocurrent around 557 pA.
  - Fig. 2B-C intensity-response: 100 ms pulses across irradiance, target K around 0.84 mW/mm2 and Imax around 642 pA.
  - Fig. 2D-F duration-response: 9.2 mW/mm2 pulses across 1-100 ms, target duration K around 3.2 ms.
  - Prolonged-current kinetics: time-to-peak around 12 ms and inactivation tau around 48 ms, but interpret this as whole-cell photocurrent kinetics rather than pure molecular ChR2 opening.
  - Fig. 3 spike output is validation only, not a photocurrent fitting target.
- Implement a NEURON soma `SEClamp` readout before changing ChR2 molecular kinetics:
  - Keep literature-style ChR2 kinetics fixed.
  - Apply ChR2 conductance/current at retained soma and dendrite segments.
  - Place `SEClamp` at soma.
  - Holding voltage: -70 mV.
  - Series/access resistance sweep: 5, 7.5, and 10 Mohm.
  - Refit `gbar` for each cell and each Rs to match 1 s / 9.2 mW/mm2 peak photocurrent around 0.557 nA.
  - Re-run the Wang protocol suite for each fitted condition.
  - Compare against the old summed-current readout to estimate how much cable/space-clamp filtering explains the fast activation mismatch.
- Current decision after first SEClamp pass:
  - Adopt soma `SEClamp` as the Wang calibration readout.
  - Use Rs = 5 Mohm as the default representative condition.
  - Keep soma depth = 100 um from the light-entry surface.
  - Do not tune light rise time yet, because both time-to-peak and intensity K are mismatched; a shutter/rise-time explanation would mainly target onset delay.
  - Next immediate sweep: effective irradiance scale while refitting `gbar` at each scale.
  - Next implementation-quality improvement: replace Python-updated `IClamp` current injection with a proper NEURON conductance mechanism for ChR2.

## ChR2 MOD Conductance Result

Implemented in:

```text
mod/chr2_4state.mod
scripts/run_wang_seclamp_mod_suite.py
outputs/wang2007_seclamp_mod_scale1/seclamp_mod_suite_results.csv
```

Fixed assumptions:

```text
readout: soma SEClamp current
Rs: 5 Mohm
soma depth from light-entry surface: 100 um
irradiance scale: 1.0
ChR2 implementation: NEURON nonspecific membrane conductance mechanism
```

Current result:

```text
PT5B_full:
  gbar = 0.07283 mS/cm2
  peak = 0.556 nA
  time to peak = 6.3 ms
  inactivation tau = 55.4 ms
  intensity K = 0.52 mW/mm2
  duration K = 2.08 ms

IT5B_full:
  gbar = 0.06364 mS/cm2
  peak = 0.560 nA
  time to peak = 6.9 ms
  inactivation tau = 55.5 ms
  intensity K = 0.49 mW/mm2
  duration K = 2.25 ms
```

Interpretation:

- Replacing Python-updated `IClamp` injection with a real NEURON membrane conductance did not substantially change the Wang suite readouts.
- The remaining mismatch is therefore not mainly caused by the Python current-injection approximation.
- With irradiance scale fixed at 1.0, the model is still too light-sensitive:
  - Wang intensity K target: about 0.84 mW/mm2
  - model: about 0.49-0.52 mW/mm2
- The model is also still too fast:
  - Wang time to peak target: about 12 ms
  - model: about 6-7 ms
  - Wang duration K target: about 3.2 ms
  - model: about 2.1-2.3 ms

Next likely model issue:

The mismatch is now more likely to be in the ChR2 light-driven transition sensitivity/rates, the effective optical calibration, or missing experimental filtering not captured by SEClamp/cable dynamics. Since the MOD implementation did not move the result much, changing the numerical implementation alone is not sufficient.

## Soma Depth Sweep Result

Implemented in:

```text
scripts/calibrate_wang2007_gbar.py
outputs/wang2007_gbar_calibration_depth_sweep/soma_entry_depth_sweep.csv
outputs/wang2007_gbar_calibration_depth_sweep/soma_entry_depth_sweep.png
```

The sweep placed the soma 50, 75, 100, 150, 200, and 250 um from the light-entry surface of a 300 um slice. The practical patch-clamp target range is treated as roughly 50-100 um below the slice surface, because very superficial cells are more likely to be cut or unhealthy, while deeper cells are harder to visualize and access.

Summary:

```text
depth from light-entry surface:
  50 um:
    PT5B gbar = 0.04715 mS/cm2, mean irradiance = 8.20 mW/mm2
    IT5B gbar = 0.04525 mS/cm2, mean irradiance = 8.06 mW/mm2
  75 um:
    PT5B gbar = 0.04342 mS/cm2, mean irradiance = 7.88 mW/mm2
    IT5B gbar = 0.04178 mS/cm2, mean irradiance = 7.75 mW/mm2
  100 um:
    PT5B gbar = 0.04228 mS/cm2, mean irradiance = 7.51 mW/mm2
    IT5B gbar = 0.03929 mS/cm2, mean irradiance = 7.45 mW/mm2
  150 um:
    PT5B gbar = 0.04165 mS/cm2, mean irradiance = 6.79 mW/mm2
    IT5B gbar = 0.03762 mS/cm2, mean irradiance = 6.79 mW/mm2
  200 um:
    PT5B gbar = 0.04191 mS/cm2, mean irradiance = 6.11 mW/mm2
    IT5B gbar = 0.03751 mS/cm2, mean irradiance = 6.13 mW/mm2
  250 um:
    PT5B gbar = 0.04429 mS/cm2, mean irradiance = 5.53 mW/mm2
    IT5B gbar = 0.04298 mS/cm2, mean irradiance = 5.62 mW/mm2
```

Recommended representative position for the next calibration step:

```text
soma depth from light-entry surface = 75-100 um
default representative value = 100 um
```

Reasoning:

- 50 um is plausible but close to the cut surface, so it may over-represent superficial, potentially damaged cells.
- 75-100 um is a common practical compromise for visualized whole-cell patching in acute slices.
- 100 um gives PT/IT estimates that are still close to the centered-slice estimate while avoiding an overly superficial soma.
- The calibrated `gbar` varies by about 10-20% over the plausible 50-150 um range, so soma depth matters but does not dominate the current-scale calibration.

## Wang Protocol Suite Result

Implemented in:

```text
scripts/calibrate_wang2007_gbar.py
outputs/wang2007_protocol_suite/wang_protocol_suite_results.csv
outputs/wang2007_protocol_suite/PT5B_full_wang_protocol_suite.png
outputs/wang2007_protocol_suite/IT5B_full_wang_protocol_suite.png
```

The current model first fits `gbar` to the Wang Fig. 2A/text peak photocurrent target:

```text
1 s pulse, 9.2 mW/mm2, target peak = 0.557 nA
```

It then tests the same `gbar` against Fig. 2B-C intensity-response and Fig. 2D-F duration-response targets.

Current results:

```text
PT5B_full:
  gbar = 0.04158 mS/cm2
  1 s / 9.2 mW/mm2 peak = 0.557 nA
  time to peak = 2.05 ms
  inactivation tau = 54.51 ms
  100 ms intensity-response K ~= 0.56 mW/mm2
  duration-response K ~= 1.50 ms

IT5B_full:
  gbar = 0.03760 mS/cm2
  1 s / 9.2 mW/mm2 peak = 0.557 nA
  time to peak = 2.10 ms
  inactivation tau = 54.56 ms
  100 ms intensity-response K ~= 0.57 mW/mm2
  duration-response K ~= 1.50 ms
```

Comparison to Wang:

```text
peak current target: 0.557 nA
  matched by construction

intensity K target: 0.84 +/- 0.2 mW/mm2
  current model: about 0.56-0.57 mW/mm2
  interpretation: too sensitive to low irradiance by about 30-35%

duration K target: 3.2 ms
  current model: about 1.5 ms
  interpretation: short-pulse activation is too fast

inactivation tau target: 48 ms
  current model: about 54.5 ms
  interpretation: reasonably close

time-to-peak target: 12 ms
  current model: about 2.1 ms
  interpretation: whole-cell current onset is too fast if the Wang 12 ms value is used literally
```

Model implication:

The current model has the right order of current scale after `gbar` fitting and a reasonable inactivation time constant. The main mismatch is early activation: it reaches peak too quickly, has too-low intensity half-activation, and has too-short duration half-activation. This points toward light-driven opening rates and/or missing experimental/cable filtering, not toward a simple expression-scale error.

## SEClamp Irradiance Scale Sweep

Implemented in:

```text
scripts/run_wang_seclamp_suite.py
outputs/wang2007_seclamp_irradiance_scale/seclamp_suite_results.csv
outputs/wang2007_seclamp_irradiance_scale/seclamp_irradiance_scale_summary.png
```

Fixed assumptions:

```text
readout: soma SEClamp current
Rs: 5 Mohm
soma depth from light-entry surface: 100 um
light rise time: ideal step, not tuned
gbar: refit at every irradiance scale to match 1 s / 9.2 mW/mm2 peak current around 0.557 nA
```

Irradiance scale results:

```text
PT5B_full:
  scale 0.4: intensity K = 1.12 mW/mm2, duration K = 2.04 ms, time to peak = 6.8 ms
  scale 0.5: intensity K = 0.91 mW/mm2, duration K = 1.99 ms, time to peak = 6.5 ms
  scale 0.6: intensity K = 0.79 mW/mm2, duration K = 1.99 ms, time to peak = 6.4 ms
  scale 0.7: intensity K = 0.71 mW/mm2, duration K = 2.03 ms, time to peak = 6.2 ms
  scale 0.8: intensity K = 0.64 mW/mm2, duration K = 2.10 ms, time to peak = 6.2 ms
  scale 1.0: intensity K = 0.52 mW/mm2, duration K = 2.22 ms, time to peak = 6.0 ms

IT5B_full:
  scale 0.4: intensity K = 0.99 mW/mm2, duration K = 2.15 ms, time to peak = 7.4 ms
  scale 0.5: intensity K = 0.83 mW/mm2, duration K = 2.11 ms, time to peak = 7.1 ms
  scale 0.6: intensity K = 0.73 mW/mm2, duration K = 2.13 ms, time to peak = 7.0 ms
  scale 0.7: intensity K = 0.65 mW/mm2, duration K = 2.19 ms, time to peak = 6.9 ms
  scale 0.8: intensity K = 0.58 mW/mm2, duration K = 2.26 ms, time to peak = 6.8 ms
  scale 1.0: intensity K = 0.49 mW/mm2, duration K = 2.39 ms, time to peak = 6.6 ms
```

Interpretation:

- Effective irradiance scaling can explain the intensity K mismatch.
- A scale around 0.5 brings PT5B and IT5B close to Wang's intensity K target of 0.84 mW/mm2.
- Effective irradiance scaling does not solve the time-to-peak mismatch. The onset remains around 6-7 ms, still faster than Wang's 12 ms.
- Effective irradiance scaling also does not solve the duration K mismatch. Duration K remains around 2.0-2.4 ms, below Wang's 3.2 ms.
- Therefore, irradiance scaling is useful but not sufficient. The remaining mismatch likely requires a more faithful conductance implementation and/or better handling of cable/space-clamp dynamics before changing molecular ChR2 rates.
