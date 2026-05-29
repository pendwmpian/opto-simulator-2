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

## Generic Q10 Wang Suite

Rationale:

- Return to Wang et al. 2007 and ignore the Chater 2010-specific fitting problem for now.
- Use a generic temperature correction on the Foutz/Nikolic-style ChR2 kinetic rates.
- Keep the interpretation intentionally modest: this is not a ChR2-specific fitted temperature model, but a sensitivity analysis using standard biological Q10 scaling.

General Q10 basis:

```text
Q10 describes how much a rate changes per 10 C temperature difference.
Ion-channel gating/inactivation rates are commonly in the approximate range Q10 = 2-4.
Classical Hodgkin-Huxley temperature scaling often uses Q10 around 3 for Na/K gating rates.
```

Sources:

```text
Temperature-robust neural function from activity dependent ion channel regulation:
  ion channel gating/inactivation Q10 values are typically around 2-4.

Kinetic Map of Kv channels:
  Hodgkin-Huxley used Q10 = 3 for Na and K conductances.
```

Experiment:

```text
baseline:
  Foutz/Nikolic-style ChR2 kinetic rates treated as 37 C reference

target:
  Wang acute-slice room-temperature condition, set h.celsius = 22 C

temperature scaling:
  rate_22 = rate_37 / Q10^((37 - 22) / 10)

Q10 sweep:
  1.5, 2.0, 2.5, 3.0

fixed assumptions:
  Rs = 5 Mohm
  soma depth = 100 um
  irradiance scale = 1.0
  ChR2 conductance implemented as NEURON .mod mechanism
  gbar refit for each Q10 to match 1 s / 9.2 mW/mm2 peak current around 0.557 nA
```

Readouts:

```text
time to peak
inactivation tau
intensity K
duration K
gbar
```

Result:

First short PT5B-only result:

```text
output:
  outputs/wang2007_mod_generic_q10_pt/seclamp_mod_suite_results.csv

cells:
  PT5B only

binary iterations:
  7, coarse first-pass only
```

```text
Q10 = 1.0:
  q10_rate_scale = 1.00
  peak = 0.591 nA
  time to peak = 4.8 ms
  inactivation tau = 55.7 ms
  intensity K = 0.41 mW/mm2
  duration K = 1.86 ms

Q10 = 2.0:
  q10_rate_scale = 0.354
  peak = 0.598 nA
  time to peak = 6.7 ms
  inactivation tau = 161.5 ms
  intensity K = 0.42 mW/mm2
  duration K = 1.91 ms

Q10 = 3.0:
  q10_rate_scale = 0.192
  peak = 0.594 nA
  time to peak = 11.1 ms
  inactivation tau = 285.2 ms
  intensity K = 0.40 mW/mm2
  duration K = 2.34 ms
```

Interpretation:

- Generic Q10 slowing can move time-to-peak toward Wang's 12 ms target.
- Q10 = 3 gets time-to-peak close, but it makes inactivation far too slow.
- Generic Q10 barely improves intensity K; it remains much lower than Wang's 0.84 mW/mm2.
- Duration K improves slightly at Q10 = 3 but remains below Wang's 3.2 ms.
- A single generic Q10 applied to all ChR2 kinetic rates is therefore not sufficient.
- This result is useful as a sensitivity check, not as a final model.

Important caveat:

This short run used only PT5B and a coarse `gbar` binary search. A production run should use more binary iterations and both PT5B/IT5B after the model choice is settled.

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

## Williams Q10 Wang Suite

Rationale:

- The generic one-Q10 experiment is too crude.
- Williams et al. 2013 / its NEURON implementation provides separate Q10 factors for H134R kinetic terms.
- We will not adopt Williams H134R absolute kinetics as the wild-type line18 model.
- We will use the Williams Q10 factors as a pragmatic temperature-scaling prior on the Foutz/Nikolic-style wild-type ChR2 scaffold.

Extracted Q10 values:

```text
Q10_Gd1       = 1.97
Q10_Gd2       = 1.77
Q10_Gr        = 2.56
Q10_e12dark   = 1.10
Q10_e21dark   = 1.95
Q10_epsilon1  = 1.46
Q10_epsilon2  = 2.77
```

Source:

```text
Williams et al. 2013 ChR2-H134R model and NEURON implementation.
The NEURON implementation applies these Q10 values to Gd1, Gd2, Gr,
e12dark, e21dark, epsilon1, and epsilon2.
```

Mapping into the current Foutz/Nikolic-style scaffold:

```text
k1 -> epsilon1 Q10
k2 -> epsilon2 Q10
gd1 -> Gd1 Q10
gd2 -> Gd2 Q10
e12 -> e12dark Q10
e21 -> e21dark Q10
gr -> Gr Q10
```

Temperature conversion:

```text
rate_22 = rate_37 / Q10^((37 - 22) / 10)
```

Experiment:

```text
Wang suite
SEClamp
Rs = 5 Mohm
soma depth = 100 um
irradiance scale = 1.0
h.celsius = 22
gbar refit per cell to match 1 s / 9.2 mW/mm2 peak current around 0.557 nA
```

Readouts:

```text
time to peak
inactivation tau
intensity K
duration K
gbar
```

Result:

First result:

```text
output:
  outputs/wang2007_mod_williams_q10/seclamp_mod_suite_results.csv

binary iterations:
  8, coarse first-pass
```

```text
PT5B_full:
  gbar = 0.08429 mS/cm2
  peak = 0.574 nA
  time to peak = 5.50 ms
  inactivation tau = 86.60 ms
  intensity K = 0.370 mW/mm2
  duration K = 1.78 ms

IT5B_full:
  gbar = 0.07041 mS/cm2
  peak = 0.561 nA
  time to peak = 6.10 ms
  inactivation tau = 86.55 ms
  intensity K = 0.369 mW/mm2
  duration K = 1.81 ms
```

Comparison to Wang:

```text
Wang peak target:
  about 0.557 nA
  model matched approximately by gbar fitting

Wang time to peak:
  about 12 ms
  model about 5.5-6.1 ms
  still too fast

Wang inactivation tau:
  about 48 ms
  model about 86.6 ms
  too slow

Wang intensity K:
  about 0.84 mW/mm2
  model about 0.37 mW/mm2
  worse than before; still too light-sensitive

Wang duration K:
  about 3.2 ms
  model about 1.8 ms
  still too fast
```

Interpretation:

- Williams-style individual Q10 scaling does not solve the Wang mismatch.
- It slows inactivation too much and makes intensity K worse.
- This reinforces that the main issue is not just temperature scaling of kinetic rates.
- The likely weak point remains the irradiance-to-transition-rate conversion / light activation formulation.

## Photon Flux Activation Formulation

Rationale:

- The previous scaffold used `ka = k * irradiance`, which directly maps mW/mm2 to transition rate.
- Williams/Foutz-style mechanisms instead convert irradiance to photon flux and then to a ChR2 transition rate using retinal cross-section and loss terms.
- This is the next required correction before interpreting intensity K.

Implementation target:

```text
Ephoton = 1e9 * h * c / wavelength_nm
flux = 1000 * irradiance_mW_mm2 / Ephoton
F = flux * sigma_retinal / (wloss * 1000)

S0 = 0.5 * (1 + tanh(120 * (100 * irradiance - 0.1)))
dp/dt = (S0 - p) / tauChR2

ka1 = epsilon1 * F * p
ka2 = epsilon2 * F * p
```

Parameters from Williams implementation:

```text
wavelength = 470 nm
hc = 1.986446e-25 J m
wloss = 1.3
sigma_retinal = 12e-20 m2
tauChR2 = 1.3 ms
epsilon1 = 0.8535
epsilon2 = 0.14
```

The absolute values come from H134R Williams implementation, so this remains a pragmatic formulation transfer. The key correction is that irradiance is now converted through photon flux instead of using arbitrary linear `k1/k2` units.

Result:

First result:

```text
output:
  outputs/wang2007_mod_photon_williams_q10/seclamp_mod_suite_results.csv

implementation caveat:
  photon flux conversion was added, but p was treated as a quasi-steady light gate.
  tauChR2 dynamic filtering is not yet active in this first pass.
```

```text
PT5B_full:
  gbar = 0.08429 mS/cm2
  peak = 0.557 nA
  time to peak = 6.90 ms
  inactivation tau = 91.97 ms
  intensity K = 0.611 mW/mm2
  duration K = 1.96 ms

IT5B_full:
  gbar = 0.07041 mS/cm2
  peak = 0.545 nA
  time to peak = 7.50 ms
  inactivation tau = 91.96 ms
  intensity K = 0.579 mW/mm2
  duration K = 2.00 ms
```

Comparison to prior Williams-Q10 linear scaffold:

```text
intensity K improved:
  linear scaffold: about 0.37 mW/mm2
  photon flux first pass: about 0.58-0.61 mW/mm2

time to peak improved modestly:
  linear scaffold: about 5.5-6.1 ms
  photon flux first pass: about 6.9-7.5 ms

inactivation tau worsened:
  linear scaffold: about 86.6 ms
  photon flux first pass: about 92 ms
```

Interpretation:

- The photon-flux conversion moves intensity K in the right direction, supporting the concern that direct `k * irradiance` was too crude.
- The model still does not match Wang:
  - Wang intensity K about 0.84 mW/mm2, model about 0.58-0.61;
  - Wang time to peak about 12 ms, model about 7 ms;
  - Wang tau about 48 ms, model about 92 ms;
  - Wang duration K about 3.2 ms, model about 2 ms.
- Remaining missing piece in this implementation: dynamic `p` / `tauChR2` filtering was simplified away to resolve the first NMODL state-equation issue. That should be restored cleanly before final interpretation.

Dynamic `p(t)` result:

```text
output:
  outputs/wang2007_mod_photon_pdynamic_williams_q10/seclamp_mod_suite_results.csv
```

The `p` gate was restored as a fifth ODE state:

```text
dp/dt = (S0 - p) / tauChR2
tauChR2 = 1.3 ms
```

Results:

```text
PT5B_full:
  gbar = 0.08429 mS/cm2
  peak = 0.554 nA
  time to peak = 8.10 ms
  inactivation tau = 91.96 ms
  intensity K = 0.607 mW/mm2
  duration K = 2.30 ms

IT5B_full:
  gbar = 0.07041 mS/cm2
  peak = 0.544 nA
  time to peak = 8.70 ms
  inactivation tau = 91.95 ms
  intensity K = 0.575 mW/mm2
  duration K = 2.33 ms
```

Interpretation:

- Restoring `p(t)` moves timing in the right direction:
  - time to peak increased from about 6.9-7.5 ms to about 8.1-8.7 ms;
  - duration K increased from about 2.0 ms to about 2.3 ms.
- Intensity K changed little and remains below Wang's 0.84 mW/mm2.
- Inactivation tau remains too slow at about 92 ms versus Wang's about 48 ms.
- Therefore, the light-gate dynamics matter, but do not fully explain Wang.

## Foutz/Nikolic Baseline Rate Correction

Problem:

The earlier `chr2_4state.mod` baseline rates were a mixed exploratory scaffold rather than a faithful Foutz/Nikolic parameter set.

Correction:

```text
gd1 = 0.130 /ms
gd2 = 0.025 /ms
e12_light = 0.053 /ms
e21_light = 0.023 /ms
e12_dark = 0.022 /ms
e21_dark = 0.011 /ms
gr = 0.0004 /ms
gamma = 0.05
```

Experiment:

```text
activation mode: photon flux
p(t): dynamic, tauChR2 = 1.3 ms
Q10 mode: Williams individual Q10 values
SEClamp Rs: 5 Mohm
soma depth: 100 um
irradiance scale: 1.0
```

Result:

```text
output:
  outputs/wang2007_mod_foutz_rates_photon_pdynamic_williams_q10/seclamp_mod_suite_results.csv

PT5B_full:
  gbar = 0.10090 mS/cm2
  peak = 0.570 nA
  time to peak = 6.50 ms
  inactivation tau = 22.07 ms
  intensity K = 0.987 mW/mm2
  duration K = 2.15 ms

IT5B_full:
  gbar = 0.08429 mS/cm2
  peak = 0.554 nA
  time to peak = 6.90 ms
  inactivation tau = 22.07 ms
  intensity K = 0.923 mW/mm2
  duration K = 2.18 ms
```

Interpretation:

- Correcting the Foutz/Nikolic baseline rates substantially improved intensity K:
  - Wang target: about 0.84 mW/mm2
  - model: about 0.92-0.99 mW/mm2
- It also fixed the previous "tau too slow" problem, but overshot:
  - Wang tau: about 48 ms
  - model: about 22 ms
- Time to peak and duration K remain too fast:
  - Wang time to peak: about 12 ms, model about 6.5-6.9 ms
  - Wang duration K: about 3.2 ms, model about 2.15-2.18 ms
- This is a much more coherent baseline than the previous mixed scaffold.

## ChR2 Localization Sweep Plan

Purpose:

Estimate how strongly Wang-suite readouts depend on assumed ChR2 membrane localization before moving to in vivo 37 C patterned stimulation.

Baseline model to use:

```text
ChR2 model:
  Foutz/Nikolic baseline rates
  photon-flux activation
  dynamic p(t)
  Williams individual Q10 scaling for 22 C

Wang suite:
  SEClamp
  Rs = 5 Mohm
  soma depth = 100 um
  irradiance scale = 1.0
  h.celsius = 22
  gbar refit per localization condition to match Wang peak photocurrent
```

Localization conditions:

```text
uniform:
  all retained soma + dendrite membrane has weight 1.0
  existing baseline

soma-only:
  soma weight = 1.0
  all dendrite weights = 0.0

dendrite-only:
  soma weight = 0.0
  all dendrite weights = 1.0

soma + proximal dendrite enriched:
  soma weight = 1.0
  proximal dendrite weight = 1.0
  distal dendrite weight = 0.25

soma + proximal dendrite reduced:
  soma weight = 0.25
  proximal dendrite weight = 0.25
  distal dendrite weight = 1.0
```

Initial proximal/distal definition:

```text
distance from soma midpoint <= 150 um:
  proximal dendrite

distance from soma midpoint > 150 um:
  distal dendrite
```

This threshold is intentionally simple for the first pass. If localization has a large effect, later sweeps can test 100/200/300 um cutoffs and apical-vs-basal localization.

Implementation plan:

1. Add a `localization_weight` to each retained segment row.
2. In `run_wang_seclamp_mod_suite.py`, multiply segment `gbar` by:

```text
effective_gbar = fitted_global_gbar * localization_weight
```

3. Keep all optical and kinetic parameters fixed.
4. For each localization condition, refit global `gbar` to Wang peak photocurrent.
5. Run the same Wang protocol suite.
6. Save one summary CSV and one comparison plot.

Readouts:

```text
fitted global gbar
retained weighted membrane area
peak current
time to peak
inactivation tau
intensity K
duration K
```

Additional diagnostic if feasible:

```text
soma current contribution
proximal dendrite contribution
distal dendrite contribution
```

Expected interpretation:

- If localization strongly changes time-to-peak and duration K, it is a meaningful parameter for later in vivo simulations.
- If localization mainly changes fitted `gbar` but not kinetics/readouts, then uniform localization is probably sufficient for the next phase.
- Soma/proximal enrichment is expected to reduce filtering delay and may make responses faster.
- Dendrite-only or distal-enriched expression is expected to require larger `gbar` and may slow apparent clamp current.

Result:

First PT5B-only result:

```text
output:
  outputs/wang2007_localization_sweep_pt/seclamp_mod_suite_results.csv
  outputs/wang2007_localization_sweep_pt/localization_summary.png

caveat:
  binary iterations = 7, so peak matching is coarse.
  This is a first-pass effect-size check, not a final calibrated run.
```

```text
uniform:
  weighted area = 19910 um2
  gbar = 0.106 mS/cm2
  peak = 0.594 nA
  time to peak = 6.5 ms
  tau = 22.1 ms
  intensity K = 0.98 mW/mm2
  duration K = 2.14 ms

soma-only:
  weighted area = 1553 um2
  gbar = 0.914 mS/cm2
  peak = 0.576 nA
  time to peak = 6.0 ms
  tau = 21.5 ms
  intensity K = 1.17 mW/mm2
  duration K = 2.15 ms

dendrite-only:
  weighted area = 8884 um2
  gbar = 0.533 mS/cm2
  peak = 0.548 nA
  time to peak = 4.3 ms
  tau = 23.6 ms
  intensity K = 0.78 mW/mm2
  duration K = 1.85 ms

soma + proximal dendrite enriched:
  weighted area = 6762 um2
  gbar = 0.311 mS/cm2
  peak = 0.529 nA
  time to peak = 6.1 ms
  tau = 22.5 ms
  intensity K = 0.96 mW/mm2
  duration K = 2.11 ms

soma + proximal dendrite reduced:
  weighted area = 6284 um2
  gbar = 0.914 mS/cm2
  peak = 0.505 nA
  time to peak = 4.2 ms
  tau = 23.1 ms
  intensity K = 0.78 mW/mm2
  duration K = 1.81 ms
```

Interpretation:

- Localization strongly changes fitted `gbar` because weighted membrane area changes.
- Localization also changes timing and intensity K, so it is not equivalent to arbitrary `gbar` scaling.
- Soma-only and proximal-enriched conditions move intensity K higher.
- Dendrite-only and proximal-reduced conditions move intensity K lower and make apparent responses faster in this current setup.
- This faster dendrite-only result is counterintuitive if one expects distal dendritic filtering to dominate; likely explanation is that the illuminated/retained dendrite population contains many relatively proximal, high-irradiance dendritic compartments and excludes soma area that contributes slower clamp load.
- A refined localization diagnostic should split soma/proximal/distal current contributions explicitly before overinterpreting the direction.

## Alternative gbar calibration: saturating photocurrent target

Rationale:

The 1 s / 9.2 mW/mm2 Wang peak photocurrent value is broad:

```text
557 +/- 153 pA
```

The intensity-response fit reports a narrower saturating current:

```text
Imax = 642 +/- 38 pA
```

Because the exact lamp/filter/objective irradiance calibration is not fully recoverable from the paper, a useful control is to fit `gbar` to the saturated photocurrent instead of the 9.2 mW/mm2 prolonged-pulse peak. This does not solve the light-delivery uncertainty, but it asks whether the expression scale inferred from a saturated response is consistent with the 9.2 mW/mm2 response.

Implementation:

```text
script:
  scripts/run_wang_seclamp_mod_suite.py

new option:
  --calibration-target imax_saturating

calibration condition:
  irradiance = 100 mW/mm2
  pulse duration = 100 ms
  target current = 0.642 nA

validation condition:
  9.2 mW/mm2
  1 s pulse
  same fitted gbar

intensity-response K calculation:
  use Wang-like 7 intensity points only:
    0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2 mW/mm2
  fit a Hill saturation curve:
    I = Imax * x^n / (K^n + x^n)
  do not include the artificial 100 mW/mm2 calibration point in this fit

other settings:
  22 degC
  Rs = 5 Mohm
  Williams Q10 mapping
  photon-flux activation
  uniform ChR2 localization
  soma depth = 100 um
```

Output:

```text
outputs/wang2007_imax_saturating_calibration/seclamp_mod_suite_results.csv
outputs/wang2007_imax_saturating_calibration/seclamp_mod_suite_results.json
```

Result:

```text
PT5B:
  fitted gbar = 0.10291 mS/cm2
  fitted saturating current = 0.642 nA
  predicted 9.2 mW/mm2 / 1 s peak = 0.580 nA
  time to peak = 6.5 ms
  inactivation tau = 22.1 ms
  intensity K = 1.08 mW/mm2
  intensity Imax from 7-point Hill fit = 0.646 nA
  intensity Hill n = 1.00
  duration K = 2.15 ms

IT5B:
  fitted gbar = 0.08942 mS/cm2
  fitted saturating current = 0.641 nA
  predicted 9.2 mW/mm2 / 1 s peak = 0.585 nA
  time to peak = 6.9 ms
  inactivation tau = 22.1 ms
  intensity K = 1.01 mW/mm2
  intensity Imax from 7-point Hill fit = 0.646 nA
  intensity Hill n = 1.00
  duration K = 2.18 ms
```

Interpretation:

- Fitting to `Imax = 642 pA` only modestly increases `gbar` relative to the previous 9.2 mW/mm2 peak-current fit.
- The same `gbar` predicts a 9.2 mW/mm2 peak current around 0.58 nA, which is close to Wang's 0.557 nA mean and well inside the reported 9.2 mW/mm2 spread.
- This supports using the saturating-current fit as the more stable expression-scale anchor.
- The intensity K should be interpreted from the Wang-like 7-point Hill fit only. The artificial 100 mW/mm2 calibration point is used to set `gbar`, not to estimate K.
- The early-kinetics mismatch remains: time-to-peak is still around 6.5-6.9 ms, duration K remains around 2.15-2.18 ms, and inactivation tau remains around 22 ms.
- Therefore the current-scale problem and the kinetic-shape problem should be separated. Saturating-current calibration can anchor `gbar`, but it does not by itself explain the faster model kinetics.
