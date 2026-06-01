# Leak Conductance Input-Resistance Sweep Plan

## Goal

Test whether reducing passive leak conductance can move the slice-pruned PT5B model toward Wang et al. 2007 Table 1 intrinsic properties.

This experiment is ChR2-free. It is a baseline-cell-state calibration experiment, not an optogenetic calibration.

## Biological Motivation

Wang et al. recorded line 18 cortical L5 pyramidal neurons across P13-P36. Younger cortical pyramidal neurons often have higher input resistance than mature neurons because of smaller effective membrane area and lower total resting conductance.

The current Dura-Bernal PT5B model appears mature and has low input resistance:

```text
full morphology Rin:
  ~56 MOhm

slice-pruned shallow soma, soma_entry_depth = 25 um:
  ~75 MOhm

Wang target:
  ChR2+ Rin = 126 +/- 13 MOhm
  ChR2- Rin = 110 +/- 4.5 MOhm
```

Because PT5B `e_pas` is strongly hyperpolarized, reducing `g_pas` should both increase input resistance and shift rest Vm in the depolarizing direction. This is biologically interpretable as a younger / lower-leak / reduced visible membrane conductance state.

## Model

```text
cell: PT5B_full
biophysics: original Dura-Bernal / NetPyNE mechanisms
ChR2: absent
temperature: 22 degC
slice geometry: vertical 300 um slice
soma_entry_depth: 25 um by default
cut face: sealed, no shunt
```

## Sweep

Primary sweep:

```text
g_pas_scale:
  1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2
```

`g_pas_scale` multiplies the original `pas.g` in every retained compartment. It does not alter active conductance density.

## Measurement

Measure input resistance with a small somatic current step:

```text
test current: 0.05 nA
step duration: 300 ms
Rin = steady delta Vm / delta I
```

Also run a short f-I check at selected currents:

```text
0, 0.1, 0.2, 0.3, 0.4, 0.5 nA
```

This checks whether `200 pA` becomes capable of producing Wang-like spiking after Rin/rest are moved toward Table 1.

## Outputs

```text
outputs/leak_conductance_input_resistance/
  leak_conductance_input_resistance.csv
  leak_conductance_input_resistance.png
```

## Interpretation

A useful calibration should move toward:

```text
Rin:
  110-126 MOhm

rest Vm:
  -58 to -61 mV
```

without producing spontaneous spikes at zero injected current or pathological depolarization block.

If reducing `g_pas` increases Rin but rest Vm remains too hyperpolarized, then passive leak alone is insufficient and HCN/hd or another depolarizing resting current needs to be considered separately.

