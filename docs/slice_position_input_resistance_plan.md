# Slice Position Input-Resistance Experiment Plan

## Goal

Evaluate whether acute-slice geometry can explain the mismatch between Wang et al. 2007 Table 1 input resistance and the current PT5B model.

The target reference is Wang line 18 cortical layer V pyramidal neurons:

```text
ChR2+ input resistance: 126 +/- 13 MOhm
ChR2- input resistance: 110 +/- 4.5 MOhm
ChR2+ resting potential: -58.3 +/- 1.8 mV
```

This experiment does not include ChR2. It only tests the electrical effect of removing membrane area by slice geometry.

## Question

If the soma is closer to the slice surface, more dendrite is cut away. Because cut dendrite should seal rather than remain open, the model should remove the lost membrane area without adding a leak at the cut face.

The key question is:

```text
Can slice-position-dependent membrane-area loss move PT5B input resistance toward ~126 MOhm?
```

## Model Assumptions

```text
cell: PT5B_full
biophysics: original Dura-Bernal / NetPyNE mechanisms
temperature: 22 degC
ChR2: absent
synaptic input: absent
slice thickness: 300 um
slice axis: z by default
cut face: sealed end, no extra shunt
```

The slice pruning should use the same geometric convention as the Wang voltage-clamp calibration:

```text
soma_entry_depth_um:
  distance from light-entry / slice surface to soma

slice_center:
  soma_axis + slice_thickness / 2 - soma_entry_depth
```

Segments outside the slice are electrically removed by setting their membrane mechanisms and capacitance near zero. This is an approximation to dendrite loss with sealed cut ends.

## Measurement

Use a small somatic current step without ChR2:

```text
baseline: no holding current
test current: 0.05 or 0.1 nA
step duration: 300 ms for fast Rin scans
temperature: 22 degC
```

For each slice position, compute:

```text
rest_v_mV
steady_v_mV during current step
apparent_input_resistance_MOhm = delta_v_mV / delta_i_nA
retained_area_um2
retained_area_fraction
spike_count during step
```

If the test current evokes a spike, that condition should be excluded from linear Rin interpretation and re-run with a smaller current.

## Sweep

First sweep:

```text
soma_entry_depth_um:
  25, 50, 75, 100, 125, 150, 200, 250

slice_axis:
  z

truncate_radius_um:
  500
```

Optional follow-up:

```text
slice_axis:
  x and z

truncate_radius_um:
  none, 300, 500
```

## Interpretation

Expected direction:

```text
smaller soma_entry_depth_um
  -> soma closer to surface
  -> more dendrite removed
  -> smaller membrane area
  -> larger input resistance
```

If Rin approaches Wang's 110-126 MOhm range by plausible slice positions, then the current-clamp Fig. 3 simulation should use the corresponding slice-pruned morphology.

If Rin remains too low even after strong pruning, the remaining mismatch likely comes from active resting conductances such as HCN/hd, passive leak density, or spike-initiation/channel-balance differences.

## Output

```text
outputs/slice_position_input_resistance/
  slice_position_input_resistance.csv
  slice_position_input_resistance.png
```

