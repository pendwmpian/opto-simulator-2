# Phase 3 Optical Model Explanation

## What Was Added For Separated Pixels

The original examples used the vertical `x-y` slice at:

```text
z = 0 um
```

For the separated-pixel pattern, many active pixels are at approximately:

```text
z = -50 um and z = +50 um
```

So additional slice outputs were generated:

```text
outputs/phase3_irradiance_zminus50/separated_pattern_xy_slice.png
outputs/phase3_irradiance_zplus50/separated_pattern_xy_slice.png
```

These slices pass through the active separated pixels and show the direct high-irradiance columns more clearly.

## What Parameters Are Currently Modeled?

The current analytical model includes:

1. Surface irradiance

```text
I0 = 1 mW/mm2
```

2. Depth attenuation

```text
exp(-mu_eff * depth)
mu_eff = 4 mm^-1
```

3. Lateral spread in the cortical surface plane

```text
sigma(depth) = 8 um + 0.18 * depth_um
```

4. Surface pixel size

```text
20 x 20 um
```

The model does not explicitly represent:

- refractive index mismatch
- phase
- polarization
- coherent interference
- wavelength-dependent Mie scattering from individual cells
- separate absorption coefficient `mu_a`
- separate scattering coefficient `mu_s`
- anisotropy factor `g`
- layer-by-layer optical heterogeneity
- blood vessel absorption

So this is not "Rayleigh scattering." Rayleigh scattering applies to particles much smaller than the wavelength. Brain tissue scattering at blue/green wavelengths is dominated by heterogeneous structures with sizes comparable to or larger than the wavelength, so Mie-like/anisotropic scattering and empirical tissue optical coefficients are more relevant.

## What Calculation Is Being Done?

The pipeline uses a Green's-function style approximation:

```text
surface pattern = sum of active square pixels
each pixel -> depth-attenuated Gaussian point-spread function
total irradiance = sum of all pixel contributions
```

For each active pixel centered at `(x0, z0)`:

```text
I(x, z, depth) =
    I0
    * exp(-mu_eff * depth)
    * pixel_area / (2*pi*sigma(depth)^2)
    * exp(-((x-x0)^2 + (z-z0)^2) / (2*sigma(depth)^2))
```

The Gaussian is normalized so that, before depth attenuation, the integrated power from each active pixel is approximately conserved while spreading laterally.

This is computationally cheap:

```text
O(number_of_active_pixels * number_of_sample_points)
```

It is therefore appropriate for rapid pattern exploration and later for computing `local_irradiance` at NEURON segment midpoints.

## Relation To Literature

The depth attenuation follows the common effective attenuation formulation described in optogenetics optics reviews:

```text
I = I0 * exp(-mu_eff * t)
mu_eff = sqrt(3 * mu_a * (mu_a + mu_s'))
mu_s' = mu_s * (1 - g)
```

where:

- `mu_a`: absorption coefficient
- `mu_s`: scattering coefficient
- `g`: anisotropy factor
- `mu_s'`: reduced scattering coefficient

Relevant references:

1. Yizhar et al., 2011, Neuron, "Optogenetics in Neural Systems"
   - Summarizes blue-light attenuation in brain tissue and optogenetic light delivery constraints.

2. Aravanis et al., 2007, Journal of Neural Engineering, "An optical neural interface..."
   - Measured blue-light transmission through rodent cortex and used a Kubelka-Munk style model.

3. Stujenske, Spellman, and Gordon, 2015, Cell Reports
   - Used Monte Carlo modeling for realistic in vivo optogenetic light and heat propagation.

4. "Optical Techniques in Optogenetics", 2015
   - Explicitly gives the effective attenuation coefficient formulation above.

5. "Scattering of Sculpted Light in Intact Brain Tissue", Scientific Reports, 2015
   - Directly relevant to this project because it uses a 488 nm laser and spatial light modulation/sculpted light.
   - Uses Monte Carlo modeling and Lorenz-Mie theory approximations to study how sculpted light degrades in intact brain tissue.

6. "Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications", eNeuro, 2016
   - Compares analytical and numerical models and discusses limitations of homogeneous brain-tissue assumptions.

## How Is This Different From Monte Carlo?

Monte Carlo light transport simulates many photon packets.

Each photon packet repeatedly undergoes:

```text
sample step length from mu_a + mu_s
move photon
absorb some weight
scatter into a new direction sampled from a phase function, often Henyey-Greenstein with anisotropy g
repeat until photon exits or weight is negligible
```

Monte Carlo can represent:

- source geometry
- angular emission distribution
- multiple scattering
- anisotropy
- absorption
- tissue boundaries
- heterogeneous optical properties
- backscatter and out-of-plane spread
- non-Gaussian long tails

The current analytical model collapses all of that into:

```text
one depth attenuation term + one Gaussian lateral blur
```

The practical differences are:

- The analytical model is much faster and easy to fit.
- It likely underestimates complex tails, backscatter, and layer/blood effects.
- It forces lateral spread to be Gaussian even if the real profile is asymmetric or heavy-tailed.
- It cannot predict heating or photon path-length distributions.
- It is sufficient for early pattern-to-segment mapping, but not sufficient for quantitative optical dosimetry.

## Recommendation

Use the current analytical model for Phase 3 and initial ChR2 coupling.

When experimental calibration data or a high-confidence tissue optical parameter set is available, replace or validate it with either:

1. measured point-spread functions through cortex at 488 nm, or
2. Monte Carlo simulation using `mu_a`, `mu_s`, `g`, refractive index, and source geometry.

