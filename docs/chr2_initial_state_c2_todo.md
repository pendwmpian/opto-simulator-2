# ChR2 Initial-State / C2 Adaptation Sweep

## Motivation

The current Wang-suite model shows a small fast current immediately after light onset. This component is mainly the fast `C1 -> O1` opening component in the 4-state ChR2 model.

Wang et al. 2007 traces do not show a clear separate fast onset component. One plausible reason is that the experimental cell was not fully dark-adapted before each pulse. Prior blue illumination, fluorescence observation, or previous pulses could leave part of ChR2 in the desensitized closed state `C2`.

The key idea is:

```text
fully dark-adapted:
  C1 = 1.0
  C2 = O1 = O2 = 0.0

partly light-adapted:
  C1 = 1 - C2_initial
  C2 = C2_initial
  O1 = O2 = 0.0
```

This avoids adding dark open conductance. It is therefore more compatible with Wang's observation that ChR2-expressing neurons did not show large baseline electrical abnormalities.

## Implementation

The NEURON mechanism now exposes initial state parameters:

```text
c1_init
o1_init
o2_init
c2_init
```

The first-pass sweep uses:

```text
script:
  scripts/sweep_chr2_initial_c2.py

cell:
  PT5B_full

temperature:
  22 degC

optical model:
  current in vitro slice model
  vitro_mu_eff = 2.12 mm^-1
  soma depth = 100 um

ChR2 kinetics:
  Foutz/Nikolic baseline rates
  Williams Q10 mapping
  photon-flux activation

gbar:
  fit once at C2_initial = 0 using saturating current target
  target = 0.642 nA at 100 mW/mm2, 100 ms

C2 sweep:
  0, 0.05, 0.10, 0.20, 0.30, 0.50

intensity protocol:
  Wang-like 7 points:
    0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2 mW/mm2
```

## Primary Readout

The most important readout is whether the small immediate light-onset response is suppressed.

For this first pass:

```text
early peak window:
  0-2 ms after light onset

early peak fraction:
  early peak current / total peak current
```

Additional readouts:

```text
9.2 mW/mm2 peak current
9.2 mW/mm2 time-to-peak
9.2 mW/mm2 inactivation tau
Wang 7-point intensity K
Hill n
```

## First-Pass Result

Output:

```text
outputs/wang2007_c2_initial_sweep/c2_initial_sweep.csv
outputs/wang2007_c2_initial_sweep/c2_initial_sweep.png
outputs/wang2007_c2_initial_sweep/c2_initial_intensity_k.png
```

Result:

```text
C2_initial = 0.00:
  9.2 mW/mm2 peak = 0.580 nA
  time-to-peak = 6.5 ms
  tau = 45.5 ms
  early peak fraction = 0.265
  intensity K = 1.08 mW/mm2
  Hill n = 1.00

C2_initial = 0.05:
  9.2 mW/mm2 peak = 0.554 nA
  time-to-peak = 6.6 ms
  tau = 46.0 ms
  early peak fraction = 0.261
  intensity K = 1.09 mW/mm2
  Hill n = 1.00

C2_initial = 0.10:
  9.2 mW/mm2 peak = 0.527 nA
  time-to-peak = 6.6 ms
  tau = 46.6 ms
  early peak fraction = 0.256
  intensity K = 1.11 mW/mm2
  Hill n = 1.00

C2_initial = 0.20:
  9.2 mW/mm2 peak = 0.473 nA
  time-to-peak = 6.8 ms
  tau = 48.0 ms
  early peak fraction = 0.245
  intensity K = 1.14 mW/mm2
  Hill n = 1.01

C2_initial = 0.30:
  9.2 mW/mm2 peak = 0.417 nA
  time-to-peak = 6.9 ms
  tau = 49.5 ms
  early peak fraction = 0.233
  intensity K = 1.19 mW/mm2
  Hill n = 1.01

C2_initial = 0.50:
  9.2 mW/mm2 peak = 0.300 nA
  time-to-peak = 7.3 ms
  tau = 54.1 ms
  early peak fraction = 0.197
  intensity K = 1.55 mW/mm2
  Hill n = 0.89
```

## Interpretation

- Increasing `C2_initial` suppresses the immediate onset component, but only modestly unless the C2 fraction is very large.
- The fast component is reduced from about 26.5% of the peak at `C2_initial = 0` to about 23.3% at `C2_initial = 0.3`, and about 19.7% at `C2_initial = 0.5`.
- Increasing `C2_initial` also reduces peak current because fewer channels begin in the readily activatable `C1` state.
- With fixed `gbar`, large `C2_initial` therefore worsens current-scale matching.
- Time-to-peak changes only weakly, from about 6.5 ms to 7.3 ms across the full sweep.
- Inactivation tau moves toward Wang's reported 48 ms around `C2_initial = 0.2-0.3`.
- Intensity K increases with C2 fraction, moving away from Wang's 0.84 mW/mm2 target.

## Current Conclusion

`C2_initial` is a plausible contributor to suppressing the early onset component, but it is not sufficient by itself unless a very large fraction of channels is pre-adapted. A very large `C2_initial` also reduces peak current and worsens intensity K.

The most useful next test is not to treat `C2_initial` as a free fit parameter by itself. Instead, use a protocol-based initial state:

```text
simulate pre-pulses or weak background light
allow a defined dark recovery interval
use the resulting C1/C2/O1/O2 state as the next trial's initial condition
```

This would connect the initial state to an experimental history rather than assigning an arbitrary C2 fraction.

## Next TODO

1. Re-run the C2 sweep with `gbar` refit at each C2 value to isolate waveform-shape effects from current-scale loss.
2. Add a protocol-adapted initial-state simulation:
   - one or more 100 ms blue pulses;
   - dark recovery intervals of 1, 2, 5, 10, and 30 s;
   - optional weak background blue-light exposure.
3. Compare three cases:
   - fully dark-adapted;
   - arbitrary C2 fraction;
   - protocol-adapted state.
4. Keep `O1_initial` and `O2_initial` at zero unless there is evidence for dark current, because open-state initialization would imply baseline conductance.
