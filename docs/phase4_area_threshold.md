# Phase 4: Centered Area Threshold

## Goal

Estimate the minimum centered stimulation area that can trigger a spike in the PT5B full morphology model.

This is the first spatial response map step. The stimulus is centered over the soma in the surface x-z plane. Location is fixed; only area changes.

## Current Assumptions

Cell:

```text
PT5B_full_cellParams.pkl
```

Optics:

```text
Yona/near-488 nm defaults
surface irradiance = 1 mW/mm2
mu_eff = 2.12 mm^-1
sigma(depth) = 8 um + 0.10 * depth_um
source sampling pitch = 2 um
```

Stimulation:

```text
488 nm equivalent irradiance
20 ms pulse
centered square or centered circle
size search range = 2-200 um
binary search iterations = 8
```

ChR2:

```text
expression scale = 1x
g_ChR2 density = 0.08 mS/cm2 per local 1 mW/mm2
E_ChR2 = 0 mV
current approximation = IClamp at segment midpoint
```

Membrane model:

```text
full PT5B morphology and topology
simplified HH/passive membrane dynamics
no inhibitory input
no synaptic background
```

## Results

Centered area thresholds:

```text
square side length threshold ~= 85.53 um
circle diameter threshold     ~= 97.13 um
```

Near-threshold ChR2 photocurrent:

```text
square 85.53 um:
  peak ChR2 current ~= 0.1884 nA
  ChR2 charge over 20 ms ~= 0.00377 nC
  first spike latency ~= 24.13 ms from pulse onset

circle 97.13 um:
  peak ChR2 current ~= 0.1889 nA
  ChR2 charge over 20 ms ~= 0.00378 nC
  first spike latency ~= 25.10 ms from pulse onset
```

At 200 um:

```text
square:
  peak current ~= 0.468 nA
  latency ~= 11.9 ms

circle:
  peak current ~= 0.435 nA
  latency ~= 12.25 ms
```

The square threshold is smaller than the circle diameter threshold because a square of side `s` has larger area than a circle of diameter `s`.

## Generated Outputs

```text
outputs/phase4_area_threshold/area_threshold_results.json
outputs/phase4_area_threshold/area_threshold_results.csv
outputs/phase4_area_threshold/area_threshold_summary.png
outputs/phase4_area_threshold/area_threshold_traces.png
```

## Interpretation

Under the current optical and ChR2 assumptions, a single 2 um spot is far too weak to drive spiking. A centered illuminated region on the order of 90-100 um is needed at 1 mW/mm2 and 1x expression.

The near-threshold square and circle conditions converge to similar total ChR2 photocurrent, about 0.19 nA. This suggests that, for centered stimuli, total recruited dendritic photocurrent is a useful first summary variable.

## Next Step

Keep the area near or above threshold, then move the stimulus center across the surface plane to build a spatial response map:

```text
stimulus center offset -> segment irradiance -> ChR2 current -> soma Vm / spike latency
```

