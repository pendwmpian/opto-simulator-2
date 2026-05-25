# Phase 3 Irradiance Pipeline

## Goal

Compute a 3D irradiance field from a surface light pattern before applying any ChR2 or membrane model.

The current pipeline answers:

```text
given a surface pattern, what irradiance reaches each x-y-z point after scattering and attenuation?
```

This is intentionally separated from ChR2 expression and cellular response. The output of this phase will later become `local_irradiance(segment, t)` for the ChR2 model.

## Initial Optical Model

The first implementation uses a simple analytical point-spread approximation.

For each ON surface pixel:

```text
I(x, z, depth) =
    I0
    * exp(-mu_eff * depth)
    * pixel_area / (2*pi*sigma(depth)^2)
    * exp(-((x-x0)^2 + (z-z0)^2) / (2*sigma(depth)^2))
```

where:

```text
I0 = 1 mW/mm2
mu_eff = 2.12 mm^-1
sigma(depth) = 8 um + 0.10 * depth_um
pixel size = 20 x 20 um
pattern field = 200 x 200 um
```

This is not yet a Monte Carlo tissue optics model. It is a testable approximation with:

- exponential depth attenuation
- lateral x-z blurring that increases with depth
- approximate conservation of source power before depth attenuation

## Example Patterns

Three example patterns were generated:

1. `single`: one 20 x 20 um pixel on.
2. `block`: contiguous 3 x 3 pixels on.
3. `separated`: five separated pixels on.

For each, the figure shows:

- left: top-view input pattern on the cortical surface
- red line: location of the vertical x-y slice
- right: irradiance in that x-y slice

Generated outputs:

```text
outputs/phase3_irradiance/single_pattern_xy_slice.png
outputs/phase3_irradiance/block_pattern_xy_slice.png
outputs/phase3_irradiance/separated_pattern_xy_slice.png
outputs/phase3_irradiance/*_depth_summary.json
outputs/phase3_irradiance/optical_params.json
```

## Quick Numerical Check

The first exploratory figures used an older conservative setting:

```text
mu_eff = 4 mm^-1
sigma(depth) = 8 um + 0.18 * depth_um
```

The code defaults have since been changed to the Yona/near-488 nm literature-based attenuation setting above. Regenerate the figures before using them quantitatively.

Depth summaries should show the expected behavior:

- At the surface, single-pixel peak is about 1 mW/mm2.
- At 300 um depth, the single-pixel peak is about 0.005 mW/mm2.
- At 300 um depth, the 3 x 3 block peak is about 0.042 mW/mm2 because neighboring pixels sum after lateral spread.
- Separated pixels do not sum strongly near the center until deeper layers where the lateral spread becomes broad.

## Next Step

The next implementation should sample this irradiance field at each NEURON segment midpoint:

```text
surface pattern -> 3D irradiance field -> segment local irradiance -> ChR2 conductance/current
```

Then spatial sweeps can ask whether soma, apical dendrite, basal dendrite, or mixed compartments dominate the spike latency and threshold.
