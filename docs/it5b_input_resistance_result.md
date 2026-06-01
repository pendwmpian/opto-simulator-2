# IT5B Input-Resistance and Current-Injection Result

## Purpose

Check whether Dura-Bernal `IT5B_full` is closer than `PT5B_full` to Wang et al. 2007 layer V cortical pyramidal neuron intrinsic properties.

This experiment excludes ChR2.

## Slice-Position Input Resistance

Script:

```text
scripts/run_slice_position_input_resistance.py
```

Output:

```text
outputs/slice_position_input_resistance_it5b/slice_position_input_resistance.csv
outputs/slice_position_input_resistance_it5b/slice_position_input_resistance.png
```

Model:

```text
cell: IT5B_full
temperature: 22 degC
test current: 50 pA
slice axis: z
slice thickness: 300 um
truncate radius: 500 um
```

Results:

| soma entry depth (um) | retained area (um2) | retained fraction | Rin (MOhm) | rest Vm (mV) |
|---:|---:|---:|---:|---:|
| full | n/a | n/a | 148.0 | -80.2 |
| 25 | 15289.2 | 0.548 | 237.1 | -80.4 |
| 50 | 18534.8 | 0.665 | 204.6 | -80.2 |
| 75 | 20118.2 | 0.722 | 191.5 | -80.1 |
| 100 | 21440.8 | 0.769 | 181.7 | -80.0 |
| 125 | 22118.9 | 0.793 | 177.0 | -80.0 |
| 150 | 22509.4 | 0.807 | 174.5 | -80.0 |
| 200 | 22716.2 | 0.815 | 173.2 | -80.0 |
| 250 | 19935.6 | 0.715 | 192.4 | -80.1 |

## Current-Injection Check

Script:

```text
scripts/run_pt5b_current_injection.py
```

Output:

```text
outputs/it5b_current_injection_22c_slice_pruned/current_injection_fi.csv
outputs/it5b_current_injection_22c_slice_pruned/current_injection_fi.png
```

Condition:

```text
cell: IT5B_full
soma_entry_depth: 100 um
temperature: 22 degC
```

Results:

| current (nA) | spikes/s | first spike latency (ms) |
|---:|---:|---:|
| 0.05 | 0 | n/a |
| 0.10 | 0 | n/a |
| 0.15 | 9 | 50.55 |
| 0.20 | 13 | 31.60 |
| 0.25 | 16 | 23.45 |
| 0.30 | 18 | 18.75 |
| 0.40 | 24 | 13.55 |
| 0.50 | 29 | 10.65 |

## Interpretation

`IT5B_full` is very different from `PT5B_full`.

Compared with Wang Table 1:

```text
Wang ChR2+ Rin:
  126 +/- 13 MOhm

IT5B full Rin:
  148 MOhm

IT5B slice-pruned Rin:
  173-237 MOhm
```

IT5B has higher input resistance than Wang, while PT5B had lower input resistance than Wang. This suggests Wang's recorded L5 pyramidal neurons may be electrically closer to an IT-like model than the Dura-Bernal PT5B model, at least for current-clamp excitability.

However, IT5B resting voltage remains much more hyperpolarized:

```text
Wang rest Vm:
  about -58 to -61 mV

IT5B rest Vm:
  about -80 mV
```

Despite the deep rest Vm, IT5B fires much more readily than PT5B:

```text
IT5B:
  0.15 nA -> 9 spikes/s
  0.40 nA -> 24 spikes/s

PT5B slice-pruned:
  0.30 nA -> 1 spike/s
  0.50 nA -> 22 spikes/s
```

## Practical Conclusion

If the goal is to reproduce Wang Fig. 3 current-clamp spiking, `IT5B_full` is a better candidate than `PT5B_full`.

If the goal is to model M1 L5B PT-like output neurons, `PT5B_full` remains more relevant biologically, but Wang Fig. 3 should not be used as a strict PT5B spiking target.

For ChR2 expression order, Wang Fig. 2 remains useful only after choosing which intrinsic cell state is being calibrated. The IT/PT difference is large enough that a single Wang-derived `gbar` should be treated as a line18 prior, not a definitive M1 PT5B expression estimate.

