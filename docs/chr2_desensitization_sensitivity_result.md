# ChR2 Desensitization Sensitivity With Electrode Readout

## Purpose

After adding the electrode/amplifier readout model, the immediate fast onset component is largely suppressed and the current scale/intensity K are close to Wang et al. 2007. The remaining mismatch is prolonged kinetics:

```text
time-to-peak too fast
inactivation tau too fast
duration K too low
```

This test checks whether small literature-near changes in the 4-state desensitization parameters can move the model toward Wang.

## Baseline Literature Values

The model starts from the Nikolic/Foutz/Grossman-style WT ChR2 4-state values:

```text
e12_light = 0.053 /ms
e21_light = 0.023 /ms
gd2       = 0.025 /ms
gamma     = 0.05
```

Parameters are changed as scale factors:

```text
e12   = e12_literature * e12_scale
e21   = e21_literature * e21_scale
gd2   = gd2_literature * gd2_scale
gamma = gamma_literature * gamma_scale
```

## Conditions Tested

All conditions use:

```text
PT5B_full
22 degC
Rs = 10 Mohm
Cpipette = 100 pF
electrode tau = 1 ms
4-stage amplifier low-pass, tau = 0.25 ms per stage
gbar fit to 100 mW/mm2, 100 ms, Imax = 0.642 nA
```

### Baseline

```text
e12_scale = 1.0
e21_scale = 1.0
gd2_scale = 1.0
gamma_scale = 1.0
```

Result:

```text
filtered peak at 9.2 mW/mm2 = 0.553 nA
filtered time-to-peak = 8.2 ms
filtered tau = 22.8 ms
early peak fraction = 0.036
intensity K = 0.817 mW/mm2
duration K = 2.00 ms
```

### Gamma-Only Increase

```text
gamma_scale = 2.0
```

Result:

```text
filtered peak at 9.2 mW/mm2 = 0.549 nA
filtered time-to-peak = 8.25 ms
filtered tau = 22.7 ms
early peak fraction = 0.035
intensity K = 0.797 mW/mm2
duration K = 2.02 ms
```

Interpretation:

Increasing `gamma` alone does not meaningfully slow prolonged kinetics.

### Slow-Desensitization Combination

```text
e12_scale = 0.5
e21_scale = 2.0
gd2_scale = 0.5
gamma_scale = 2.0
```

Result:

```text
filtered peak at 9.2 mW/mm2 = 0.582 nA
filtered time-to-peak = 8.65 ms
filtered tau = 29.9 ms
early peak fraction = 0.034
intensity K = 0.589 mW/mm2
duration K = 2.08 ms
```

Interpretation:

- Slowing desensitization in this way moves inactivation tau in the right direction.
- The effect is not large enough to reach Wang's about 48 ms tau.
- The same change makes intensity K too low, moving away from Wang's about 0.84 mW/mm2.
- Duration K remains around 2 ms and does not approach Wang's about 3.2 ms.

## Current Conclusion

Small literature-near changes to `e12/e21/gd2/gamma` are not enough to solve the prolonged kinetics mismatch cleanly.

The fast onset component is best explained by recording-chain filtering. The remaining prolonged kinetics mismatch may require one of:

```text
1. a different ChR2 kinetic model component for slow adaptation/desensitization;
2. voltage-dependent ChR2 kinetics not captured by the current simplified model;
3. more complete experimental protocol history and light adaptation;
4. a more realistic clamp/electrode feedback model if it affects slow current readout more than expected.
```

The current tested parameter changes should not be adopted as defaults because they improve tau only partially while degrading intensity K.
