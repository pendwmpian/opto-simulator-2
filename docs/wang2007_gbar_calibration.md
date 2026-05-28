# Wang 2007 Line 18 gbar Calibration

## Purpose

Estimate an effective uniform ChR2 `gbar` for Thy1-ChR2-YFP line 18 before collecting our own patch data.

This calibration uses Wang et al., 2007 acute-slice photocurrent measurements as the anchor, rather than choosing an arbitrary `1x` or `4x` expression scale.

## Experimental Conditions Parsed from Wang 2007

Reference:

```text
Wang et al., 2007, PNAS
High-performance genetically targetable optical neural silencing by light-driven proton pumps
and optogenetic activation using Thy1-ChR2-YFP line 18 data.
```

General slice and recording conditions:

```text
preparation: acute cortical parasagittal slices
slice thickness: 250-350 um
mouse age: P13-P36
line: Thy1-ChR2-YFP line 18
cell: ChR2-positive cortical layer V pyramidal neuron
recording: whole-cell recording with Axoclamp 2B
voltage-clamp holding voltage for photocurrent traces: -70 mV
temperature: 21-24 C
patch pipette: 2-7 Mohm
internal solution: K-gluconate based
extracellular solution: ACSF
retinal/retinol: sometimes added at 1 uM, no detectable effect
synaptic blockers for isolated photocurrents: CNQX + APV + picrotoxin or GABAzine
wide-field ChR2 excitation: 465-495 nm arc-lamp light, electronic shutter
mapping excitation: 488 or 458 nm argon laser spot
```

The paper uses several different photostimulation protocols. They should not be collapsed into a single condition.

### Fig. 2A: Example prolonged photocurrent

```text
protocol: prolonged wide-field light pulse
wavelength: 465-495 nm
area: large field, about 0.4 mm2
holding voltage: -70 mV
cell type: cortical layer V ChR2-positive pyramidal neuron
purpose: show direct ChR2 photocurrent vs wild-type negative control
reported numeric intensity in text: 9.2 mW/mm2 for 1-s pulses
figure caveat: the example trace scale bar visually suggests the displayed pulse may be below the maximum scale used elsewhere
```

The text explicitly states that 1-s wide-field pulses at 9.2 mW/mm2 produced inward currents and reports the mean peak current under those conditions. However, the visual light-monitor trace in Fig. 2A should not be used as an independent calibrated intensity unless the original figure data are available.

### Fig. 2B-C: Intensity-response photocurrents

```text
protocol: wide-field light pulses of varied intensity
duration: 100 ms
wavelength: 465-495 nm
area: large field, about 0.4 mm2
holding voltage: -70 mV
readout: peak photocurrent amplitude
fit: Hill equation
K: 0.84 +/- 0.2 mW/mm2
Imax: 642 +/- 38 pA
Hill n: 0.76 +/- 0.1
```

This is the best source for irradiance-current saturation, not for 1-s inactivation kinetics.

### Fig. 2D-F: Duration-response photocurrents

```text
protocol: wide-field light pulses of varied duration
intensity: 9.2 mW/mm2
brief duration panel: 1-8 ms
longer duration panel: 10-100 ms
readout: peak photocurrent and charge
duration K for half-maximal photoactivation: 3.2 ms
Hill n for duration-response: 7.7
```

This is the best source for how quickly photocurrent peak grows with pulse duration under strong light.

### Prolonged-illumination kinetics

```text
condition stated in text: prolonged illumination at 9.2 mW/mm2
time to peak: 12.0 +/- 0.1 ms
inactivation tau: 48 +/- 0.9 ms
n: 10
```

This target belongs to wide-field voltage-clamp photocurrent kinetics during sustained illumination. It should not be confused with action-potential latency in current clamp.

### Fig. 3A-B: Current-clamp action-potential output

```text
protocol: light flashes of varied intensity
duration for intensity-response summary: 1 s
readout: number of action potentials
mean latency to peak of first action potential: 6.1 +/- 0.1 ms
jitter: 0.3 ms
half-maximal light intensity for AP number: 0.49 +/- 0.1 mW/mm2
Hill n: 0.82 +/- 0.1
maximum AP number during 1-s flashes: 25.4 +/- 1.2
```

The 6.1 ms value is not photocurrent time-to-peak. It is current-clamp spike timing, measured as latency from light onset to the peak of the first action potential. It depends on membrane excitability, spike threshold, and current-clamp dynamics.

### Fig. 3C-D: Frequency following

```text
protocol: repeated light pulses
pulse duration: 4 ms
intensity: 9.2 mW/mm2
readout: probability of evoked action potentials
reliable following: up to about 30 Hz
50% success: about 35 Hz
```

This is useful for validating short-pulse spike reliability, but it is not a direct voltage-clamp photocurrent calibration.

### Fig. 4: Spatial sensitivity mapping

```text
protocol: scanned laser spot over single dye-filled cortical pyramidal neuron
spot diameter: 1-100 um
wavelength: 458 or 488 nm
pulse duration: 4 ms
example intensity: 48 mW/mm2
maximum intensity mentioned: up to 70 mW/mm2
readout: local depolarization / action potential map
main result: responses could be evoked over much of the cell, but APs required soma/proximal region
```

This figure is closest to the later spatial-stimulation problem, but it uses high-intensity small-spot laser stimulation and current-clamp responses, not wide-field voltage-clamp photocurrent.

## Calibration Targets

Use different targets for different calibration steps:

```text
Step A: current scale / expression
  use Fig. 2A text:
    1-s wide-field pulse
    9.2 mW/mm2
    peak photocurrent = 557 +/- 153 pA

Step B: intensity saturation
  use Fig. 2B-C:
    100-ms wide-field pulses
    K = 0.84 +/- 0.2 mW/mm2
    Imax = 642 +/- 38 pA
    Hill n = 0.76 +/- 0.1

Step C: sustained-current kinetics
  use prolonged-illumination text:
    9.2 mW/mm2
    time to peak = 12.0 +/- 0.1 ms
    inactivation tau = 48 +/- 0.9 ms

Step D: short-pulse activation
  use Fig. 2D-F:
    9.2 mW/mm2
    peak-current duration K = 3.2 ms
    duration-response Hill n = 7.7

Step E: spike-output validation, not photocurrent fitting
  use Fig. 3:
    first AP peak latency = 6.1 +/- 0.1 ms
    4-ms pulses follow to about 30 Hz
```

Do not use the Fig. 3 first-spike latency as the target for ChR2 photocurrent rise time. It is a downstream current-clamp output.

```text
mean peak photocurrent: 557 +/- 153 pA
Imax: 642 +/- 38 pA
half-max irradiance K: 0.84 +/- 0.2 mW/mm2
time to peak: 12.0 +/- 0.1 ms
inactivation tau: 48 +/- 0.9 ms
```

## Simulation Setup

Separate in vitro optical model:

```text
large-field illumination
not the Yona in vivo patterned-light model
vertical slice model
slice thickness = 300 um
slice plane = morphology x-y
slice thickness axis = morphology z
light propagation axis = morphology z
vitro illumination mode = depth_decay
vitro mu_eff = 2.12 mm^-1
truncate dendrites farther than 500 um from soma
uniform ChR2 gbar density across retained soma + dendrite membrane
axon excluded
voltage fixed at -70 mV
```

Wang et al. report parasagittal cortical slices, not coronal slices, with thickness 250-350 um. The goal here is to reproduce that acute-slice photocurrent condition, not a planned user experiment. In this morphology coordinate system, y is treated as the cortical depth axis. The model therefore uses a vertical cortical slice containing y and one horizontal axis. We choose morphology x-y as the slice face and z as the 300 um thickness axis.

Parasagittal and coronal slices are similar for this single-cell calibration in the sense that both are vertical cortical slices containing the depth axis. They are not anatomically identical: a parasagittal slice preserves dorsal-ventral and anterior-posterior axes, whereas a coronal slice preserves dorsal-ventral and medial-lateral axes. Because the Dura-Bernal morphology coordinates are not registered to a specific AP/ML brain coordinate frame, this calibration only uses the common vertical-slice geometry: depth axis plus one horizontal axis plus a finite 300 um thickness.

The slice truncation is a first approximation to acute-slice dendrite cutting. It is not a detailed slice-cut geometry. The current implementation keeps compartments inside the 300 um z-thickness slab and within 500 um of the soma.

The optical model is now slice-specific: light enters along the thickness axis and decays through the 300 um slice. The x-y irradiance figure is shown from the light-entry side of the slice. This is still approximate, because the exact Wang illumination path through the microscope objective and slice is not fully reconstructed.

## Results

Both PT5B and IT5B were calibrated to match the Wang mean peak photocurrent of 0.557 nA.

```text
PT5B_full:
  gbar = 0.04158 mS/cm2
  peak current = 0.55700 nA
  steady current = 0.27397 nA
  time to peak = 2.05 ms
  inactivation tau = 54.51 ms
  retained membrane area = 20323.5 um2
  mean irradiance on retained membrane = 7.00 mW/mm2

IT5B_full:
  gbar = 0.03760 mS/cm2
  peak current = 0.55700 nA
  steady current = 0.27418 nA
  time to peak = 2.10 ms
  inactivation tau = 54.56 ms
  retained membrane area = 22509.4 um2
  mean irradiance on retained membrane = 6.86 mW/mm2
```

The PT and IT effective `gbar` estimates are within about 17% of each other. This suggests the calibration is not extremely sensitive to the chosen L5B morphology, at least for large-field voltage-clamp current.

## Important Caveat

The fitted peak current and inactivation tau are close to the Wang target, but the simulated time-to-peak is too fast:

```text
target time to peak: 12 ms
current model: about 2 ms
```

This means the current ChR2 kinetics open too quickly. The `gbar` estimate is useful as a first order current scale, but the kinetics parameters should be tuned before using this model for precise latency predictions.

The likely next adjustment is to slow the light-driven opening rates while preserving peak current by refitting `gbar`.

Possible causes of the time-to-peak mismatch:

- The current Python kinetics are a simplified scaffold, not a direct transcription of the full Nikolic/Grossman/Foutz parameter set.
- The light-driven opening rates are probably too fast.
- The voltage-clamp simulation treats the entire retained membrane as immediately clamped to -70 mV; real dendritic voltage clamp is imperfect and can slow whole-cell current onset.
- The exact Wang illumination waveform may have finite shutter/filter rise time rather than an ideal step.
- Room-temperature kinetics and wavelength/spectrum differences can shift rates.
- The reported 12 ms may reflect whole-cell current filtering and distributed dendritic charging, not just molecular ChR2 opening.

The next calibration should fit kinetics in this order:

1. Tune opening rates to match time-to-peak around 12 ms.
2. Tune desensitization/transition rates to match inactivation tau around 48 ms.
3. Refit `gbar` to match peak current around 557 pA.
4. Keep the fitted `gbar` as `line18_wang2007_default`.

## Generated Outputs

```text
outputs/wang2007_gbar_calibration_vertical_slice/wang2007_protocol.json
outputs/wang2007_gbar_calibration_vertical_slice/calibration_results.json
outputs/wang2007_gbar_calibration_vertical_slice/calibration_results.csv
outputs/wang2007_gbar_calibration_vertical_slice/PT5B_full_vitro_irradiance_xy.png
outputs/wang2007_gbar_calibration_vertical_slice/PT5B_full_vitro_irradiance_xz.png
outputs/wang2007_gbar_calibration_vertical_slice/PT5B_full_slice_light_sideview.png
outputs/wang2007_gbar_calibration_vertical_slice/IT5B_full_vitro_irradiance_xy.png
outputs/wang2007_gbar_calibration_vertical_slice/IT5B_full_vitro_irradiance_xz.png
outputs/wang2007_gbar_calibration_vertical_slice/IT5B_full_slice_light_sideview.png
outputs/wang2007_gbar_calibration_vertical_slice/PT5B_full_calibrated_current.png
outputs/wang2007_gbar_calibration_vertical_slice/IT5B_full_calibrated_current.png
```
