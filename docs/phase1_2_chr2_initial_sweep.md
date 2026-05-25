# Phase 1-2 Initial ChR2 Sweep

## Scope

This is the first executable single-cell optogenetic sweep using the Dura-Bernal PT5B full morphology.

Important limitation: this run does not yet use the original Dura-Bernal NMODL channel set. The full `pt3d` morphology and topology are used, but membrane dynamics are currently simplified to NEURON built-in HH/passive channels. ChR2 is approximated as a light-dependent inward current injected into each illuminated segment using `IClamp`.

This makes the result useful for first-order geometry, irradiance, attenuation, and soma-vs-dendrite contribution checks. It is not yet a final biophysical ChR2 model.

## Model

Cell:

```text
external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl
```

Light:

```text
wavelength: 488 nm
surface irradiance sweep: 0.1, 0.2, 0.5, 1, 2, 5 mW/mm2
depth attenuation: exp(-depth / 250 um)
pulse: 20 ms, starting at 100 ms
```

ChR2 expression sweep:

```text
0.25x, 1x, 4x
```

Illumination targets:

```text
soma
dendrite = apical + basal
apical
basal
all = soma + apical + basal, excluding axon
```

## Approximation Used For ChR2

For each segment:

```text
local_irradiance = surface_irradiance * exp(-depth / attenuation_length)
g_ChR2 = base_g * expression_scale * local_irradiance * segment_area
I_ChR2 ~= g_ChR2 * (E_ChR2 - V_ref)
```

Current implementation:

```text
E_ChR2 = 0 mV
V_ref = -70 mV
base_g = 0.08 mS/cm2 per mW/mm2
```

Because this is injected as `IClamp`, it does not yet update with instantaneous membrane voltage. A conductance-based MOD mechanism should replace this in the next iteration.

## Current Illumination Assumption

The current simulation assumes the cell is fully inside a uniformly illuminated surface patch.

Therefore:

- `mW/mm2` is treated as surface irradiance.
- Every non-axon segment receives light according to its membrane area and depth attenuation.
- Larger illuminated dendritic area produces proportionally larger total ChR2 current.
- No finite 2 um pixel edge, lateral scattering, or spatial spot profile is included yet.

In short, this is a wide-field illumination model over the entire cell, not a single-pixel stimulation model.

## Initial Findings

The first sweep suggests that, under this attenuation and expression assumption, dendritic illumination dominates over soma-only illumination.

Key observations:

- Soma-only stimulation does not spike even at 5 mW/mm2 and 4x expression.
- At 1x expression, dendrite/apical/all stimulation begins to spike around 1 mW/mm2.
- At 4x expression, dendrite/apical/all stimulation spikes around 0.2 mW/mm2.
- Basal-only stimulation is weaker at low-to-mid irradiance but can spike at high expression and 5 mW/mm2.
- Latency shortens as irradiance increases. In the 1x dendrite condition, latency goes from about 13 ms at 1 mW/mm2 to about 6.6 ms at 5 mW/mm2.
- Apical-only and all-dendrite conditions are similar in threshold in this geometry, implying that superficial apical dendrites are currently carrying most of the light-driven effect.

These findings are plausible for surface illumination because apical dendrites are closer to the surface and receive less attenuated light than soma or basal dendrites.

With `mu_eff = 4 mm^-1` and `I0 = 1 mW/mm2`, the estimated non-axon light interception is:

```text
apical dendrite: 91.3%, mean local irradiance 0.373 mW/mm2
basal dendrite:   7.5%, mean local irradiance 0.049 mW/mm2
soma:             1.2%, mean local irradiance 0.052 mW/mm2
```

This explains why soma-only stimulation is weak while apical/dendritic stimulation can trigger spikes.

## Generated Outputs

```text
outputs/phase1_chr2/sweep_results.json
outputs/phase1_chr2/sweep_expr_0.25.png
outputs/phase1_chr2/sweep_expr_1.png
outputs/phase1_chr2/sweep_expr_4.png
outputs/phase1_chr2/example_traces.png
outputs/light_distribution/PT5B_light_mu4_I1.png
outputs/light_distribution/PT5B_light_mu4_I5.png
outputs/light_distribution/PT5B_light_summary_mu4_I1.csv
outputs/light_distribution/PT5B_light_summary_mu4_I5.csv
```

## Next Corrections

Before treating this as quantitative:

1. Compile and load the original Dura-Bernal NMODL mechanisms.
2. Instantiate the PT5B full cell with its original channel densities.
3. Replace the `IClamp` ChR2 approximation with a conductance-based ChR2 mechanism.
4. Calibrate 488 nm tissue attenuation with literature or measured power-through-tissue values.
5. Add explicit spatial spot profiles, not only depth attenuation.
6. Compare continuous illumination and PWM illumination at matched total light dose.
