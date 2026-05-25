# Optical Parameter Provenance

## Current Code Parameters

The current Phase 3 code uses:

```text
surface_irradiance_mw_mm2 = 1.0
mu_eff_mm_inv = 2.12
lateral_sigma0_um = 8.0
lateral_spread_per_depth = 0.10
pixel_size = 20 x 20 um
```

These are now the deterministic defaults. The depth attenuation is tied to near-488 nm cortical optical properties described below. The lateral spread parameter is still an analytical PSF approximation and should be replaced by a Yona-style beam-spread function or Monte Carlo model if exact optical transport becomes necessary.

- literature-based blue-light depth attenuation
- lateral blurring that increases with depth
- summation of nearby illuminated pixels
- fast enough computation for pattern exploration

## Literature Values Identified So Far

### 488 nm or near-488 nm optical properties

Wide-field optical mapping / Monte Carlo literature reports a 488 nm cortical parameter set:

```text
lambda = 488 nm
mu_a = 0.33 mm^-1
mu_s = 21 mm^-1
g = 0.8
mu_s' = mu_s * (1 - g) = 4.2 mm^-1
```

Source:

```text
Wide-field optical mapping of neural activity and brain haemodynamics:
considerations and novel approaches
```

The source states these values assume 3% blood volume in tissue and 75% average oxygen saturation.

### 473 nm mouse cortex / optogenetics measurements

Yona et al., eNeuro 2016, used mouse cortical slices illuminated from above by a 473 nm laser from a 200 um, NA 0.37 fiber and compared analytical and Monte Carlo models.

Reported / discussed values include:

```text
published brain-tissue scattering coefficient examples:
mu_s = 168.6 cm^-1 = 16.86 mm^-1 at 453 nm
mu_s = 120 cm^-1 = 12.0 mm^-1 at 480 nm

fitted cortical parameters:
mu_s = 211 cm^-1 = 21.1 mm^-1
g = 0.86
mu_s' = 21.1 * (1 - 0.86) = 2.95 mm^-1

estimated cortical scattering length:
~47 um
```

This supports strong blue-light scattering in cortex, but does not by itself give a full 488 nm `mu_a, mu_s, g` set for our exact window/patterned projection geometry.

### Sculpted 488 nm light

The Scientific Reports 2015 paper "Scattering of Sculpted Light in Intact Brain Tissue, with implications for Optogenetics" is especially relevant because it uses:

```text
488 nm laser
spatial light modulator
sculpted/focused light
Monte Carlo scattering model
Lorenz-Mie-supported single-scattering approximation
```

However, the experimental preparation is zebrafish larval brain, not mouse cortex. It is therefore useful for qualitative expectations about sculpted-light degradation, but should not be copied directly as mouse cortical optical coefficients.

## How Current `mu_eff = 2.12 mm^-1` Relates To Literature

If using the diffusion-style effective attenuation expression:

```text
mu_eff = sqrt(3 * mu_a * (mu_a + mu_s'))
```

then the 488 nm parameter set:

```text
mu_a = 0.33 mm^-1
mu_s' = 4.2 mm^-1
```

gives:

```text
mu_eff = sqrt(3 * 0.33 * (0.33 + 4.2))
       ~= 2.12 mm^-1
```

The current:

```text
mu_eff = 2.12 mm^-1
```

matches this homogeneous-tissue diffusion estimate.

Previous exploratory runs used `mu_eff = 4 mm^-1`, which was intentionally conservative. That value should not be treated as the default unless patch data cannot be explained by ChR2 expression/sensitivity after the optical model is fixed.

Recommended optical handling:

```text
fix optical parameters from literature or Monte Carlo
fit ChR2 expression/sensitivity to patch data
```

## Blood Vessel Absorption

The current model includes only uniform bulk absorption via `mu_eff`. It does not explicitly model vessels.

This is acceptable for an initial average-field model, especially if the illuminated region avoids large surface vessels. It is not acceptable for predicting local pixel-to-dendrite coupling near a visible vessel.

Expected influence:

- Large pial vessels can cast strong local shadows at 488 nm because hemoglobin absorption is substantial in the blue-green range.
- Capillaries contribute more like a spatially distributed background absorption if no individual vessel map is available.
- Hemodynamic state and oxygenation can change the effective absorption, although for optogenetic stimulation this is typically slower than the immediate ChR2 response.

Practical model options:

1. Uniform blood volume approximation:

```text
use bulk mu_a that already assumes blood volume and oxygenation
```

2. Vessel mask approximation:

```text
I_surface(x,z) <- I_surface(x,z) * exp(-mu_a_blood * vessel_path_length(x,z))
```

3. Full Monte Carlo voxel model:

```text
assign each voxel tissue/blood optical properties and simulate photon packets
```

For this project, option 1 is enough for early ChR2 parameter fitting. Option 2 becomes worthwhile if surface vessel images are available from the cranial window.

## Neural Cell Scattering

The current L5B morphology cannot by itself predict the tissue scattering field.

Reason:

- The optical field is shaped by many cells, nuclei, membranes, myelin, extracellular space, glia, blood, and tissue interfaces.
- One reconstructed L5B cell is a negligible fraction of the optical volume.
- The current morphology lacks refractive index maps, nuclei, organelles, and surrounding tissue geometry.

What can be computed with the L5B morphology:

```text
given an optical irradiance field, sample irradiance at each dendritic/somatic segment
```

What cannot be computed from the L5B morphology alone:

```text
how that cell scatters the incident light field in a realistic cortical volume
```

To model cell-body-rich versus neuropil-rich scattering, use literature or measured optical coefficients. The 488 nm sculpted-light paper explicitly notes that cell nuclei can contribute strongly to scattering, but that has to be represented as tissue optical parameters or voxel heterogeneity, not from a single electrical morphology.

## Recommended Simulation Path

Because local light intensity cannot be measured directly in the planned experiments, use patch responses for calibration. To avoid overfitting and identifiability problems, keep the optical parameters fixed to literature/Monte Carlo values by default and fit ChR2 expression/sensitivity first.

1. Keep optical model fixed by default:

```text
Yona/near-488 nm optical parameters or Monte Carlo-derived PSF
```

2. Keep ChR2 expression parameterized:

```text
g_ChR2_density / expression_scale
```

3. Fit against patch observables:

```text
subthreshold depolarization
spike threshold irradiance
first spike latency
adaptation during repeated pulses
soma vs dendrite-targeted pattern differences
```

4. Treat optical and ChR2 parameters as partially confounded:

```text
stronger ChR2 expression and stronger local irradiance can produce similar voltage responses
```

5. Only revisit optical parameters if ChR2 expression cannot explain:

```text
spatial pattern dependence
soma-biased vs apical-biased response differences
latency scaling with irradiance
depth-dependent response trends
```

6. Break the confound using spatial patterns:

```text
soma-biased vs apical-biased vs basal-biased patterns
```
