# Chater 2010 ChR2 Temperature Calibration TODO

## Purpose

Estimate a defensible 22 C ChR2 kinetics parameter set without fitting kinetics to Wang et al. 2007.

The Wang suite should remain a line18 expression/current-scale validation where `gbar` is fit. ChR2 molecular kinetics should instead be constrained by independent ChR2 biophysics data.

## Source Logic

- Use the Foutz/Nikolic wild-type ChR2 four-state model as the baseline photocycle structure.
- Treat the Foutz/Nikolic parameter set as the high-temperature/physiological simulation reference because Foutz used it in a 37 C neuronal simulation context.
- Use Chater et al. 2010 as the independent 22 C macroscopic ChR2 kinetics constraint.
- Do not fit ChR2 kinetics to Wang et al. 2007.

## Chater 2010 Targets

At approximately 22 C:

```text
2 ms flash:
  time to peak = 2.9 +/- 0.1 ms
  decay tau = 7.2 +/- 0.7 ms

1 ms flash:
  deactivation tau = 10.7 +/- 0.4 ms

300 ms stimulus:
  deactivation tau = 11.5 +/- 0.5 ms
  desensitization approximately 75%

Recovery from desensitization:
  half-time = 3.1 s at -85 mV
```

## Model Degrees Of Freedom

Do not freely fit every rate. Use grouped kinetic scales:

```text
s_open:
  light-driven opening sensitivity/rate

s_close:
  O1/O2 closing rates

s_adapt:
  O1/O2 adaptation transition rates

s_recovery:
  recovery from desensitized closed state
```

This keeps the number of adjustable dimensions close to the number of independent Chater constraints.

## Simulation Plan

1. Build a morphology-free ChR2-only simulator.
2. Use normalized current; no `gbar` fitting.
3. Simulate Chater-style protocols:
   - 1 ms flash
   - 2 ms flash
   - 300 ms stimulus
   - paired-pulse recovery after a 200 ms pre-desensitizing pulse
4. Grid-search grouped rate scales.
5. Save the best 22 C candidate parameter set.
6. Only after this, inject the chosen 22 C parameter set into the Wang SEClamp suite and fit `gbar`.

## Results

First coarse grid completed:

```text
script:
  scripts/calibrate_chr2_chater2010.py

outputs:
  outputs/chater2010_chr2_temperature/chater2010_grid_results.csv
  outputs/chater2010_chr2_temperature/chater2010_best.json
  outputs/chater2010_chr2_temperature/chater2010_best_traces.png
```

Best coarse grouped-scale result:

```text
s_open = 0.7
s_close = 1.0
s_adapt = 0.5
s_recovery = 1.0

2 ms flash time to peak = 1.98 ms
2 ms flash decay tau = 6.98 ms
1 ms flash decay tau = 7.00 ms
300 ms stimulus off-decay tau = 13.43 ms
300 ms desensitization fraction = 0.56
recovery half-time = 3.05 s
```

Comparison to Chater targets:

```text
2 ms flash time to peak:
  target = 2.9 ms
  model = 1.98 ms
  interpretation = opening / peak timing is still too fast

2 ms flash decay tau:
  target = 7.2 ms
  model = 6.98 ms
  interpretation = good

1 ms flash deactivation tau:
  target = 10.7 ms
  model = 7.00 ms
  interpretation = too fast

300 ms stimulus off-decay tau:
  target = 11.5 ms
  model = 13.43 ms
  interpretation = reasonably close

300 ms desensitization:
  target = about 75%
  model = 56%
  interpretation = insufficient desensitization

recovery half-time:
  target = 3.1 s
  model = 3.05 s
  interpretation = good
```

Interpretation:

- A simple grouped scale of the Foutz/Nikolic-style four-state scaffold can reproduce recovery and some decay kinetics.
- It does not yet reproduce all Chater 22 C constraints simultaneously.
- The main remaining issues are:
  - peak timing after a 2 ms flash is too early;
  - 1 ms flash deactivation is too fast;
  - 300 ms desensitization is too weak.
- This suggests the current scaffold still lacks some details of the published Nikolic/Foutz photon-flux formulation, especially the explicit `T`, `Lambda1/Lambda2`, and light/dark transition handling.

Next steps:

1. Port the Foutz/Nikolic equations more faithfully instead of using the older simplified `k1_per_mw_ms/k2_per_mw_ms` scaffold.
2. Then rerun the same Chater 22 C evaluation.
3. Only if the faithful model still cannot match Chater with grouped temperature scaling should we consider a different ChR2 model family.

## Williams Q10 Transfer Experiment

Rationale:

- Williams et al. 2013 models ChR2(H134R), not wild-type ChR2.
- The absolute Williams parameters should not replace the wild-type Foutz/Nikolic parameter set.
- However, the temperature scaling / Q10 values from Williams can be used as a plausible ChR2-family temperature correction.
- To avoid overclaiming, sweep Q10 around the Williams-inspired range and validate against Chater 2010 wild-type 22 C kinetics.

Experiment:

```text
baseline:
  Foutz/Nikolic wild-type ChR2 rates treated as 37 C reference

temperature conversion:
  rate_22 = rate_37 / Q10^((37 - 22) / 10)

Q10 sweep:
  1.5, 2.0, 2.5, 3.0

rate groups:
  opening rates
  closing rates
  adaptation/desensitization rates
  recovery rate

first pass:
  apply same Q10 to all kinetic groups

second pass:
  optionally allow group-specific Q10 values around the same range if the one-Q10 model fails badly
```

Validation:

```text
Use the same Chater targets:
  2 ms flash time to peak = 2.9 ms
  2 ms flash decay tau = 7.2 ms
  1 ms flash decay tau = 10.7 ms
  300 ms stimulus off-decay tau = 11.5 ms
  300 ms desensitization ~= 75%
  recovery half-time = 3.1 s
```

Decision rule:

- If one-Q10 scaling gets close enough, use that 22 C parameter set for Wang suite.
- If one-Q10 scaling fails selectively, inspect which group needs a different Q10.
- Do not fit to Wang kinetics.

First Q10 transfer result:

```text
outputs:
  outputs/chater2010_q10_transfer/chater2010_q10_results.csv
  outputs/chater2010_q10_group_transfer/chater2010_q10_group_results.csv
```

One-Q10 scaling result:

```text
Q10 values tested: 1.5, 2.0, 2.5, 3.0
best one-Q10 value: 1.5

2 ms flash time to peak = 1.98 ms
2 ms flash decay tau = 10.01 ms
1 ms flash decay tau = 10.14 ms
300 ms stimulus off-decay tau = 19.14 ms
300 ms desensitization fraction = 0.54
recovery half-time = not reached in the tested interval range
```

Interpretation:

- A single Q10 applied to all kinetic rates is not adequate.
- It improves 1 ms flash deactivation but makes longer-stimulus decay too slow and does not solve the too-early 2 ms flash peak.

Group-specific Q10 result:

```text
Q10 candidates per group: 1.0, 1.5, 2.0, 2.5
best:
  q10_open = 1.0
  q10_close = 1.0
  q10_adapt = 2.5
  q10_recovery = 1.5

2 ms flash time to peak = 1.98 ms
2 ms flash decay tau = 7.18 ms
1 ms flash decay tau = 7.20 ms
300 ms stimulus off-decay tau = 11.68 ms
300 ms desensitization fraction = 0.59
recovery half-time = 4.03 s
```

Interpretation:

- Group-specific Q10 improves 2 ms decay and 300 ms off-decay.
- It still does not reproduce Chater fully:
  - 2 ms flash peak remains too early;
  - 1 ms flash deactivation remains too fast;
  - 300 ms desensitization remains too weak;
  - recovery is slower than target.
- Best result leaves opening and closing unscaled, while strongly temperature-scaling adaptation. This is not a clean mechanistic solution.

Current conclusion:

Q10 transfer alone is not enough with the current simplified light-opening scaffold. The next scientifically cleaner step is to port the full Foutz/Nikolic photon-flux activation formula before drawing conclusions about whether Williams-style Q10 transfer can explain Chater 22 C kinetics.
