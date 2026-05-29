# Wang 2007 Fig. 2 Protocol and Protocol-History Assessment

## Purpose

Before adding more ChR2 states or fitting prolonged kinetics, separate Wang et al. 2007 Fig. 2 readouts by protocol and decide whether protocol history / light adaptation can be modeled from the paper.

## Fig. 2 Readouts

### Fig. 2A / Text: Prolonged Photocurrent

Use for:

```text
current scale under sustained wide-field illumination
prolonged whole-cell photocurrent kinetics
```

Known conditions:

```text
preparation: acute cortical parasagittal slice
cell: ChR2-positive cortical layer V pyramidal neuron
recording: whole-cell voltage clamp
holding voltage: -70 mV
light: 465-495 nm arc-lamp, electronic shutter
wide-field area: about 0.4 mm2
stimulation: 1 s pulse
reported high intensity: 9.2 mW/mm2
reported peak current: 557 +/- 153 pA
reported time-to-peak: 12.0 +/- 0.1 ms
reported inactivation tau: 48 +/- 0.9 ms
```

Caveat:

The figure's light-monitor trace and the text do not fully expose all stimulus history. The example trace should not be overused as a calibrated trace beyond the values stated in the text.

### Fig. 2B-C: Intensity-Response

Use for:

```text
irradiance-current saturation
Imax
Hill K
Hill n
```

Known conditions:

```text
duration: 100 ms
light: 465-495 nm wide-field
readout: peak photocurrent
fit: Hill equation
K = 0.84 +/- 0.2 mW/mm2
Imax = 642 +/- 38 pA
Hill n = 0.76 +/- 0.1
```

Likely intensity points:

```text
0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2 mW/mm2
```

Caveat:

The paper does not state the order of intensities, inter-pulse interval, number of repeated sweeps per cell, or whether cells were fully dark-adapted before each pulse.

### Fig. 2D-F: Duration-Response

Use for:

```text
short-pulse photocurrent recruitment
how peak current grows with pulse duration at strong irradiance
```

Known conditions:

```text
intensity: 9.2 mW/mm2
durations: 1-8 ms and 10-100 ms panels
readout: peak photocurrent and charge
duration K = 3.2 ms
duration Hill n = 7.7
```

Caveat:

The duration-response K is not a current time-to-peak. It is the pulse duration needed to recruit half-maximal peak current under strong light.

## Protocol-History / Light-Adaptation Information

Information needed to model protocol history mechanistically:

```text
order of intensities in Fig. 2B
order of durations in Fig. 2D-F
inter-pulse intervals
number of repeats and averaging
dark adaptation time before each protocol
whether fluorescence/YFP observation light was present before stimulation
whether low-intensity pulses were always preceded by high-intensity pulses
whether sweeps were randomized or monotonic
```

Information available in the paper:

```text
not enough for protocol-history reconstruction
```

Interpretation:

The paper contains enough information to simulate the stimulus waveforms themselves, but not enough information to set a defensible trial-by-trial ChR2 initial state from protocol history.

## Modeling Decision

Do not add a protocol-history or light-adapted initial-state model as a default Wang calibration step.

Reason:

```text
Any specific pre-pulse history, dark recovery interval, or background-light condition would be speculative.
```

Allowed use:

```text
sensitivity analysis only
```

Example:

```text
fully dark-adapted state:
  C1 = 1

light-adapted sensitivity:
  C2_initial sweep

protocol-history simulation:
  only if an explicit hypothetical protocol is stated
```

But these should not be used to tune the final line18 `gbar` or kinetics unless independent experimental metadata become available.

## Consequence for Next Simulations

Use Wang Fig. 2 as separate validation modules:

```text
1. Expression/current scale:
   saturating Imax = 642 +/- 38 pA

2. Intensity-response validation:
   100 ms pulses
   7-point Hill fit
   K = 0.84 +/- 0.2 mW/mm2
   Hill n = 0.76 +/- 0.1

3. Duration-response validation:
   9.2 mW/mm2
   1-100 ms pulse duration series
   duration K = 3.2 ms
   Hill n = 7.7

4. Prolonged kinetics:
   1 s pulse
   9.2 mW/mm2
   time-to-peak = 12 ms
   tau = 48 ms
   treat as lower-confidence for molecular kinetic fitting because recording chain and example/protocol details are under-specified
```

Then move to Fig. 3 current-clamp validation rather than adding speculative molecular states.

## Fig. 3 Simulation Priority

Fig. 3 should be treated as downstream validation of the calibrated photocurrent model, not as a photocurrent-fitting target.

Useful targets:

```text
1 s light flashes:
  spike count vs irradiance
  half-max light intensity for AP number
  first AP peak latency

4 ms repeated pulses:
  spike probability vs frequency
  reliable following to about 30 Hz
```

This is useful because the final in vivo problem is about whether patterned light can drive spiking/motor output, not only voltage-clamp photocurrent.

## Current Recommendation

Do not proceed to a new slow-adaptation branch yet.

Instead:

```text
1. Keep the electrode/amplifier readout model because it explains the missing fast onset component and improves intensity K.
2. Keep ChR2 4-state kinetics at literature values for now.
3. Treat protocol-history/light-adaptation as sensitivity only, not calibration.
4. Next simulate Fig. 2D-F duration response and Fig. 3 current-clamp outputs with the current best model.
5. Revisit slow adaptation only if Fig. 2D-F and Fig. 3 fail in a way that specifically requires slow ChR2 recovery/adaptation.
```
