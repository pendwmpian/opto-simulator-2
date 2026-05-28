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
