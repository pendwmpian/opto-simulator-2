# Slice Position Input-Resistance Result

## Purpose

Test whether acute-slice membrane-area loss can explain the low input resistance of the PT5B current-clamp model.

This experiment excludes ChR2 and measures only the intrinsic electrical input resistance of the original Dura-Bernal PT5B cell after slice-position-dependent electrical pruning.

Reference from Wang et al. 2007 Table 1:

```text
ChR2+ L5 pyramidal neurons:
  input resistance = 126 +/- 13 MOhm
  resting potential = -58.3 +/- 1.8 mV

ChR2- neurons:
  input resistance = 110 +/- 4.5 MOhm
```

## Method

Script:

```text
scripts/run_slice_position_input_resistance.py
```

Command:

```text
python scripts/run_slice_position_input_resistance.py \
  --out-dir outputs/slice_position_input_resistance \
  --include-full \
  --target-c 22 \
  --test-current-nA 0.05 \
  --soma-entry-depths-um 25,50,75,100,125,150,200,250
```

Model:

```text
cell: PT5B_full
biophysics: original Dura-Bernal / NetPyNE mechanisms
temperature: 22 degC
ChR2: absent
synaptic input: absent
slice thickness: 300 um
slice axis: z
truncate radius: 500 um
test current: 50 pA
```

Slice pruning follows the Wang voltage-clamp geometry. Segments outside the slice are electrically removed by setting their membrane density mechanisms and capacitance near zero. No cut-face shunt is added; this approximates sealed dendritic cut ends.

## Output

```text
outputs/slice_position_input_resistance/slice_position_input_resistance.csv
outputs/slice_position_input_resistance/slice_position_input_resistance.png
```

## Results

![slice position input resistance](/Users/pend/Documents/Code/opto-simulator-2/outputs/slice_position_input_resistance/slice_position_input_resistance.png)

| soma entry depth (um) | retained area (um2) | retained fraction | Rin (MOhm) | rest Vm (mV) |
|---:|---:|---:|---:|---:|
| full | n/a | n/a | 56.3 | -77.0 |
| 25 | 15285.7 | 0.529 | 75.3 | -75.8 |
| 50 | 17771.1 | 0.615 | 68.5 | -76.2 |
| 75 | 19339.2 | 0.669 | 63.9 | -76.3 |
| 100 | 19910.2 | 0.688 | 62.1 | -76.2 |
| 125 | 20223.5 | 0.699 | 61.5 | -76.3 |
| 150 | 20323.5 | 0.703 | 61.4 | -76.3 |
| 200 | 20323.5 | 0.703 | 61.4 | -76.3 |
| 250 | 19359.8 | 0.669 | 62.9 | -76.2 |

## Interpretation

Moving the soma toward the slice surface increases input resistance in the expected direction:

```text
full morphology:
  Rin = 56.3 MOhm

soma entry depth = 25 um:
  retained area = 52.9%
  Rin = 75.3 MOhm
```

However, even very superficial soma placement does not reach Wang's 110-126 MOhm range. Slice-position-dependent dendritic loss explains part of the discrepancy but not all of it.

The resting potential also remains around -76 mV across the sweep, far from Wang's -58.3 mV. This indicates that the mismatch is not just membrane area. Active resting conductances and passive reversal/leak balance remain important.

## Conclusion

The current best slice geometry for increasing Rin is the shallowest tested condition:

```text
soma_entry_depth_um = 25
```

But this should not be treated as fully calibrated to Wang Table 1. It improves Rin from about 56 MOhm to about 75 MOhm, still below the target.

Next likely contributors:

```text
1. HCN/hd resting conductance
2. passive leak density and e_pas distribution
3. total retained morphology still larger than the actual recorded acute-slice cell
4. spike-threshold / Na availability / AIS representation
```

Before Fig. 3 ChR2 current-clamp simulation is interpreted quantitatively, the cell's current-clamp baseline should be calibrated to Wang Table 1 or a documented reason should be given for accepting lower Rin.

