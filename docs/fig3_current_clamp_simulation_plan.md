# Wang 2007 Fig. 3 Current-Clamp Simulation Plan

## Goal

Use the calibrated in vitro ChR2 model to test whether a current-clamp / voltage-recording simulation can reproduce Wang et al. 2007 Fig. 3-like spiking responses.

This is downstream validation. It should not be used to refit ChR2 photocurrent parameters unless there is a clear failure mode.

## Starting Model

Use the current best Wang photocurrent model:

```text
cell: PT5B_full first, IT5B_full optional follow-up
ChR2 expression scale: gbar from Wang saturating-current calibration
ChR2 kinetics: literature 4-state values
temperature: 22 degC
slice optical attenuation: mu_eff = 1.3 mm^-1
soma depth: 100 um
illumination: Wang-style wide-field, 465-495 nm represented by current photon-flux model
```

Current best `gbar` from voltage-clamp/electrode-readout calibration:

```text
PT5B_full:
  gbar ~= 0.13478 mS/cm2
```

This `gbar` was calibrated to:

```text
100 mW/mm2, 100 ms, saturating current target Imax = 0.642 nA
```

and validated at:

```text
9.2 mW/mm2, 1 s:
  filtered peak current ~= 0.557 nA
```

## Important Change From Voltage Clamp

Do not reuse the voltage-clamp electrode/amplifier current readout as-is.

Reason:

```text
Fig. 2 voltage clamp:
  measured current
  SEClamp/access resistance and current readout filtering are relevant

Fig. 3 current clamp:
  measured membrane potential
  current-clamp bridge/access resistance and voltage recording filter are relevant
```

For Fig. 3, the primary simulated biological readout is:

```text
soma membrane voltage
spike count
first spike peak latency
spike probability / reliability
```

The recording readout should be separated into:

```text
biophysical Vm:
  soma.v directly from the cell model

recorded Vm:
  optional low-pass filtered soma.v
  optional bridge/access voltage artifact model
```

The first implementation should report both, but use raw soma membrane voltage for spike detection unless the voltage filter is shown to alter peak timing materially.

## Fig. 3 Targets

### Fig. 3A-B: Light-Evoked Spiking

Known targets:

```text
protocol: 1 s light flashes of varied intensity
readout: number of action potentials
half-maximal light intensity for AP number: 0.49 +/- 0.1 mW/mm2
Hill n: 0.82 +/- 0.1
first AP peak latency: 6.1 +/- 0.1 ms
jitter: 0.3 ms
maximum AP number during 1-s flashes: 25.4 +/- 1.2
```

Important:

```text
first AP peak latency is not ChR2 photocurrent time-to-peak
```

It depends on:

```text
resting Vm
spike threshold
Na/K channel state
baseline current
initial membrane state
recording mode
cell health
```

### Fig. 3C-D: Frequency Following

Known targets:

```text
protocol: repeated 4 ms light pulses
intensity: 9.2 mW/mm2
readout: evoked AP probability
reliable following: up to about 30 Hz
50% success frequency: about 35 Hz
```

## Current-Clamp Assumptions

The Fig. 3 simulation requires assumptions that were not needed for voltage clamp:

```text
initial membrane potential
baseline holding current
whether the cell is at natural rest or manually biased
initial gating states of Na/K/HCN channels
series/access resistance and bridge balance
voltage-recording filtering
spike detection threshold
```

First-pass assumptions:

```text
temperature:
  h.celsius = 22

initial Vm:
  start at the model's natural resting potential after equilibration

baseline current:
  no holding current in the first pass

synaptic input:
  none

recording mode:
  current clamp with no injected current

voltage readout:
  raw soma.v plus optional low-pass displayed trace
```

If the model does not spike at Wang-like intensities, add a second pass:

```text
small DC bias current to set baseline Vm near a plausible L5 pyramidal current-clamp resting voltage
```

This should be reported as a separate condition, not hidden inside the default model.

## Implementation Steps

1. Add a Fig. 3 current-clamp script:

```text
scripts/run_wang_fig3_current_clamp.py
```

2. Reuse existing pieces:

```text
CellFromNetPyNE
prepare_template_rows
add_localization_weights
install_chr2
williams_q10_scales
mu_eff = 1.3 mm^-1
gbar = 0.13478 mS/cm2 by default
```

3. Add current-clamp simulation:

```text
insert ChR2 conductance
set h.celsius = 22
initialize/equilibrate cell
optionally apply IClamp DC bias
apply light by changing segment chr2_4state.irr over time
record soma.v and spike times
```

4. Fig. 3A-B protocol:

```text
intensity values:
  0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2 mW/mm2

pulse duration:
  1 s

readouts:
  spike count during 1 s
  first AP peak latency
  maximum firing rate / spike count saturation
  Hill fit for spike count vs intensity
```

5. Fig. 3C-D protocol:

```text
intensity:
  9.2 mW/mm2

pulse duration:
  4 ms

frequencies:
  5, 10, 20, 30, 35, 40, 50 Hz

readouts:
  evoked spike probability per pulse
  first spike latency per pulse
  failure rate
```

6. Output:

```text
outputs/wang2007_fig3_current_clamp/
  fig3_intensity_spike_count.csv
  fig3_intensity_spike_count.png
  fig3_example_voltage_traces.png
  fig3_frequency_following.csv
  fig3_frequency_following.png
```

7. Documentation:

```text
docs/fig3_current_clamp_result.md
```

## Key Risks

### Cell Model May Not Match Wang's Recorded Cells

Dura-Bernal PT5B morphology/channels are not necessarily the exact L5 pyramidal cells Wang recorded.

This affects:

```text
resting Vm
threshold
spike latency
firing frequency
adaptation
```

Therefore Fig. 3 should be treated as qualitative validation unless current-clamp baseline excitability is first calibrated.

### Baseline Vm Is Critical

If the cell starts too hyperpolarized, it may not spike even with correct ChR2 current.

If it starts near threshold, spike latency can become much shorter and spike count can increase.

Therefore report:

```text
resting Vm before light
whether DC bias was used
```

### Voltage-Clamp Recording Filter Does Not Transfer

The electrode/amplifier current filtering used for Fig. 2 explains current traces, but Fig. 3 voltage traces need their own treatment.

Use:

```text
raw soma.v for biological spike detection
optional voltage low-pass only for displayed trace
```

### Protocol History Remains Unknown

Do not add C2 initial adaptation or pulse-history effects to Fig. 3 by default. These can be sensitivity analyses later.

## First-Pass Success Criteria

The model is promising if it shows:

```text
increasing spike count with irradiance
half-max spike count intensity on the order of 0.5 mW/mm2
first AP peak latency on the order of several ms
4 ms pulses can evoke spikes up to tens of Hz
```

Exact matching is not expected until baseline excitability and current-clamp recording conditions are constrained.
