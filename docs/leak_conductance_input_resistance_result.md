# Leak Conductance Input-Resistance Sweep Result

## Purpose

Test whether reducing passive leak conductance in a slice-pruned PT5B model can move intrinsic properties toward Wang et al. 2007 Table 1.

This experiment does not include ChR2.

## Model

```text
cell: PT5B_full
biophysics: original Dura-Bernal / NetPyNE mechanisms
temperature: 22 degC
slice geometry: vertical 300 um slice
soma_entry_depth: 25 um
retained membrane area: 15285.7 um2
retained area fraction: 0.529
cut face: sealed, no shunt
```

Sweep:

```text
g_pas_scale:
  1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2
```

Output:

```text
outputs/leak_conductance_input_resistance/leak_conductance_input_resistance.csv
outputs/leak_conductance_input_resistance/leak_conductance_fi.csv
outputs/leak_conductance_input_resistance/leak_conductance_input_resistance.png
```

## Results

![leak conductance sweep](/Users/pend/Documents/Code/opto-simulator-2/outputs/leak_conductance_input_resistance/leak_conductance_input_resistance.png)

### Input Resistance

| g_pas scale | rest Vm (mV) | Rin (MOhm) |
|---:|---:|---:|
| 1.0 | -75.54 | 73.3 |
| 0.8 | -74.91 | 76.9 |
| 0.6 | -74.20 | 81.1 |
| 0.5 | -73.82 | 83.4 |
| 0.4 | -73.41 | 85.9 |
| 0.3 | -72.98 | 88.6 |
| 0.2 | -72.52 | 91.6 |

Reducing leak conductance shifts both metrics in the expected direction:

```text
g_pas down:
  Rin increases
  rest Vm becomes less hyperpolarized
```

However, even an 80% reduction in passive leak only reaches about 92 MOhm. This remains below Wang's 110-126 MOhm range.

### f-I Check

The f-I response improves as `g_pas` decreases:

| g_pas scale | 0.2 nA spikes/s | 0.3 nA spikes/s | 0.4 nA spikes/s | 0.5 nA spikes/s |
|---:|---:|---:|---:|---:|
| 1.0 | 0 | 15 | 23 | 31 |
| 0.8 | 1 | 17 | 25 | 32 |
| 0.6 | 9 | 19 | 26 | 33 |
| 0.5 | 11 | 20 | 27 | 34 |
| 0.4 | 12 | 20 | 28 | 35 |
| 0.3 | 13 | 21 | 28 | 35 |
| 0.2 | 14 | 22 | 29 | 36 |

This moves the model toward Wang Fig. 3, but `200 pA` still does not produce `20-25 spikes/s`.

## Interpretation

Reducing passive leak is biologically defensible as a young / lower-leak / lower visible membrane conductance state, especially because Wang's line 18 recordings span P13-P36. It helps, but it is not sufficient by itself.

The remaining mismatch likely requires one or more of:

```text
1. stronger effective dendrite/axon loss than the current 25 um slice-depth pruning
2. changes in HCN/hd or other resting active conductances
3. altered spike initiation / Na availability
4. a different L5 pyramidal subtype than Dura-Bernal PT5B
5. a depolarized vitro cell state not captured by passive leak reduction alone
```

## Practical Conclusion

For a Wang-vitro PT5B-like baseline, `g_pas_scale = 0.4-0.6` is a reasonable first sensitivity range:

```text
Rin:
  81-86 MOhm

rest Vm:
  -74 to -73 mV

200 pA response:
  9-12 spikes/s
```

This is still not a full Wang Table 1 match. It should be treated as a partial correction, not a final calibrated baseline.

