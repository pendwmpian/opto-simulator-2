# 488 nm Light Scattering Assumptions

## Current Illumination Geometry

The current single-cell sweep assumes wide-field surface illumination.

This means:

- The whole cell is inside the illuminated field.
- Surface irradiance is spatially uniform across the cell's x-z footprint.
- The only optical loss currently modeled is depth-dependent attenuation.
- Segment input scales with membrane area, so a larger illuminated dendritic membrane area receives proportionally more total ChR2 current.
- Axon is excluded from optical stimulation.

The current model does not yet simulate a focused 2 um pixel, finite spot edge, diffraction, speckle, or lateral redistribution from scattering. Those should be added before interpreting patterned-pixel stimulation quantitatively.

## Depth Attenuation Model

For the first light-distribution visualization, local irradiance is:

```text
I(depth) = I0 * exp(-mu_eff * depth)
```

where:

```text
I0: surface irradiance at 488 nm, in mW/mm2
depth: distance from the local apical-most point toward deeper tissue, in mm
mu_eff: effective attenuation coefficient, in mm^-1
```

The initial value used here is:

```text
mu_eff = 2.12 mm^-1
attenuation length = 1 / mu_eff = 471.7 um
```

This is a simple effective model for blue light in cortex based on near-488 nm literature values. It is still an analytical approximation; a Monte Carlo or beam-spread-function model should replace it if exact optical transport is required.

## Literature Basis

The exact value for mouse cortex at 488 nm depends on preparation, cortical depth, blood, numerical aperture, beam shape, and whether one reports ballistic, scattered, or total fluence.

The following literature supports using strong attenuation for blue light and gives the modeling framework:

1. Yizhar et al., 2011, Neuron, "Optogenetics in Neural Systems"
   - Summarizes measured light transmission in brain tissue at 473, 561, 594, and 635 nm.
   - Shows that blue light attenuates strongly with distance in brain tissue.
   - Refers to Aravanis et al. measurements and an optogenetics light power density calculator.

2. Stujenske, Spellman, and Gordon, 2015, Cell Reports, "Modeling the Spatiotemporal Dynamics of Light and Heat Propagation for In Vivo Optogenetics"
   - Provides realistic light and heat propagation modeling for in vivo optogenetic stimulation.
   - Useful for the next step when estimating tissue heating and comparing wavelength-dependent spread.

3. "Optical Techniques in Optogenetics", 2015
   - Gives the standard effective attenuation formulation:

```text
I = I0 * exp(-mu_eff * t)
mu_eff = sqrt(3 * mu_a * (mu_a + mu_s'))
mu_s' = mu_s * (1 - g)
```

4. Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications, eNeuro, 2016
   - Uses the classic 473 nm, 200 um core, NA 0.37 fiber configuration from Aravanis/Yizhar-style measurements.
   - Provides a bridge between empirical measurements and analytical/Monte Carlo models.

## Interpretation For This Project

For the current window-based 488 nm patterned stimulation, the exact source geometry differs from a fiber. Therefore, using a fiber-derived light cone directly would be inappropriate.

The current simplification is:

```text
wide-field surface pattern -> local membrane irradiance by depth attenuation
```

This matches the question "if a cell is under an illuminated patch, which compartments receive how much light?" but not yet "how much does a 2 um surface pixel spread laterally before reaching the dendrite?"

The next optical refinement should separate:

1. surface pattern blur/spread in x-z
2. depth attenuation in y
3. local membrane area integration
4. wavelength-specific calibration at 488 nm

## 1x ChR2 Expression Definition

In the current code, `1x expression` is not yet a measured Thy1-ChR2 expression density. It is a reference conductance scale:

```text
g_ChR2_density = 0.08 mS/cm2 per (mW/mm2 of local irradiance)
```

At local irradiance of 1 mW/mm2:

```text
g_ChR2_density = 0.08 mS/cm2
E_ChR2 = 0 mV
V_ref = -70 mV
current_density ~= 5.6 uA/cm2
```

For 1 um2 membrane area:

```text
area = 1e-8 cm2
I_ChR2 ~= 0.056 pA per um2 at 1 mW/mm2
```

Thus:

```text
0.25x = 0.02 mS/cm2 per mW/mm2
4x    = 0.32 mS/cm2 per mW/mm2
```

This parameter should later be fitted to patch-clamp data or literature photocurrent measurements in Thy1-ChR2 L5B cells. Until then it should be interpreted as a sensitivity scale, not a biological expression rate.
