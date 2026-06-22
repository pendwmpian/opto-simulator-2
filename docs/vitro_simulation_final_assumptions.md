# Final In Vitro Simulation Assumptions and Provenance

## Scope

This document records the assumptions used in the final Wang et al. 2007 in vitro voltage-clamp calibration branch.

It distinguishes among:

- conditions reported directly by Wang et al.;
- literature-based model choices;
- project-specific approximations;
- exploratory assumptions that were tested but not retained.

The final reference is the soma-centered `0.4 mm2` field sweep in:

```text
outputs/wang_4state_field_position_soma_fine_sweep/
```

The representative center condition produced approximately:

```text
gbar = 0.1009 mS/cm2
9.2 mW/mm2, 1 s peak current = 0.578 nA
steady/peak = 0.228
time to peak = 9.9 ms
decay tau = 27.4 ms
intensity K = 0.858 mW/mm2
duration K = 1.68 ms
```

These values should always be interpreted together with the assumptions below.

## 1. Experimental Target

### Source experiment

The calibration target is Wang et al. 2007, Thy1-ChR2-YFP line 18.

Reported experimental conditions:

| Item | Value | Status |
|---|---|---|
| Preparation | acute cortical parasagittal slice | reported by Wang |
| Slice thickness | 250-350 um | reported by Wang |
| Mouse age | P13-P36 | reported by Wang |
| Recorded cells | ChR2-positive cortical layer V pyramidal neurons | reported by Wang |
| Brain region | described as cerebral cortex; not identified as M1 PT5B | important mismatch |
| Temperature | 21-24 C | reported by Wang |
| Voltage clamp | whole-cell, holding potential -70 mV | reported by Wang |
| Amplifier | Axoclamp 2B | reported by Wang |
| Pipette resistance | 2-7 MOhm | reported by Wang |
| Wide-field source | arc lamp, band-pass 465-495 nm | reported by Wang |
| Pulse control | Uniblitz T132 electronic shutter | reported by Wang |
| Large illuminated area | approximately 0.4 mm2 | reported by Wang |
| Synaptic blockers | CNQX, APV, and picrotoxin or GABAzine were sometimes used | reported by Wang |

Wang photocurrent targets used in this project:

| Readout | Target |
|---|---:|
| Peak current at 9.2 mW/mm2, 1 s | 557 +/- 153 pA |
| Saturating intensity-response current, Imax | 642 +/- 38 pA |
| Intensity half-activation K | 0.84 +/- 0.2 mW/mm2 |
| Intensity Hill coefficient | 0.76 +/- 0.1 |
| Prolonged-stimulus time to peak | 12.0 +/- 0.1 ms |
| Inactivation tau | 48 +/- 0.9 ms |
| Duration half-activation K | 3.2 ms |
| Duration Hill coefficient | 7.7 |

The seven intensity points used for the simulated Hill fit are inferred from Fig. 2B/C as:

```text
0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2 mW/mm2
```

Their exact numeric list is not printed in the Wang Methods. It is a figure-based reconstruction.

Reference:

- Wang H, et al. High-speed mapping of synaptic connectivity using photostimulation in Channelrhodopsin-2 transgenic mice. PNAS. 2007. https://doi.org/10.1073/pnas.0700384104

## 2. Cell Model

### Adopted cell

The electrical cell is the Dura-Bernal et al. `PT5B_full` multicompartment M1 model, including its original active membrane mechanisms and conductance distributions.

Reference:

- Dura-Bernal S, et al. Multiscale model of primary motor cortex circuits predicts in vivo cell-type-specific, behavioral state-dependent dynamics. Cell Reports. 2023. https://doi.org/10.1016/j.celrep.2023.112574

### Important biological mismatch

Wang recorded generic cortical layer V pyramidal neurons from line 18 and did not establish that the Fig. 2 cells were M1 corticospinal/PT5B neurons. The Dura-Bernal PT5B cell is therefore a morphology/biophysics proxy, not a cell-identity match.

Consequences:

- voltage-clamp photocurrent amplitude remains useful for estimating an effective ChR2 conductance density;
- Wang Table 1 resting voltage, input resistance, and Fig. 3 firing should not be treated as direct PT5B validation targets;
- the final `gbar` is conditional on this PT5B morphology and conductance distribution.

### Axon

Axonal segments do not receive ChR2 in the optical rows. Any axonal structure present in the original cell may still participate electrically. No attempt is made to reconstruct the long in vivo axonal arbor lost in slice preparation.

## 3. Slice Geometry

### Adopted geometry

| Item | Model value | Basis |
|---|---:|---|
| Slice thickness | 300 um | midpoint of Wang's 250-350 um range |
| Slice type | vertical cortical slice abstraction | Wang used parasagittal slices |
| Slice face | morphology x-y plane | coordinate convention |
| Thickness/light axis | morphology z | coordinate convention |
| Soma depth from light-entry face | 100 um | plausible patching depth assumption |

The Dura-Bernal coordinates are not registered to the anatomical AP/ML axes. Therefore the model preserves the concept of a vertical slice containing the cortical-depth axis, but it does not claim exact parasagittal anatomical registration.

### No radial truncation

The earlier `truncate_radius_um = 500` cutoff was erroneous because it removed valid distal apical morphology. It is not part of the final model.

Final retained non-axonal membrane area in the slice-selection rows:

```text
approximately 26,459 um2
```

### Critical implementation limitation

The final Wang voltage-clamp sweep does not fully remove slice-excluded dendrites from the electrical cell.

Current behavior:

```text
full Dura-Bernal electrical morphology is instantiated
only segments selected inside the 300 um slab receive ChR2 conductance
segments outside the optical slice rows remain electrically present but have no ChR2
```

Thus the final model is best described as:

```text
full electrical PT5B cell + slice-restricted ChR2 illumination
```

It is not a fully cut and sealed acute-slice morphology. Separate input-resistance experiments did electrically remove out-of-slice membrane and assumed sealed cut ends, but that pruning was not propagated into the final Wang voltage-clamp sweep.

## 4. Illumination Geometry

### Field size and shape

Wang reports a large illuminated area of approximately `0.4 mm2`, but not its exact aperture shape or center.

Project assumptions:

```text
shape = circle
area = 0.4 mm2
radius = sqrt(0.4 mm2 / pi) = 356.8 um
center = patched soma
```

The circular shape and soma-centered placement are not stated by Wang. They are experimental-practice assumptions.

The soma-centered field was swept by `-100, 0, +100 um` in both slice-face axes. The resulting `gbar`, time-to-peak, intensity K, and duration K were nearly invariant. Therefore the final result is locally robust to reasonable centering error around the soma.

### Spectrum

Wang used filtered arc-lamp light spanning `465-495 nm`.

The ChR2 photon conversion uses one effective wavelength:

```text
wavelength_nm = 470
```

This ignores the arc-lamp spectral shape, filter transmission, objective transmission, and wavelength-dependent ChR2 action spectrum. Their net effect is absorbed into the effective irradiance and fitted `gbar`.

### Surface irradiance

Reported irradiance values are treated as irradiance at the tissue/light-entry surface. Water, objective, shutter, and optical-path losses are not modeled separately.

## 5. In Vitro Light Attenuation

### Adopted model

The final in vitro branch uses homogeneous exponential depth attenuation:

```text
I(z) = I_surface * exp(-mu_eff * z)
mu_eff = 1.3 mm^-1
```

Only depth attenuation is applied inside the circular field. There is no Monte Carlo photon transport and no explicit lateral scattering across the field edge.

### Provenance and interpretation

Yona et al. measured and modeled blue-light propagation in cortical tissue using analytical and Monte Carlo approaches. Their work supports strong scattering, anisotropy, and the need to distinguish tissue optical transport from a simple ray model.

Reference:

- Yona G, Meitav N, Kahn I, Shoham S. Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications. eNeuro. 2016. https://doi.org/10.1523/ENEURO.0059-15.2015

However, `mu_eff = 1.3 mm^-1` is not a direct verbatim parameter reported by Wang and should not be described as an exact Yona measurement for this preparation. It is a project-level effective attenuation chosen to represent acute slice tissue with less absorption than the earlier in vivo-oriented `mu_eff = 2.12 mm^-1`.

The earlier `2.12 mm^-1` value came from the diffusion-style relation:

```text
mu_eff = sqrt(3 * mu_a * (mu_a + mu_s'))
```

using near-488 nm bulk-tissue assumptions. It was not retained for the final in vitro branch because it may include vascular/blood absorption absent from an acute slice.

### Effects omitted

- explicit blood vessels and blood-flow-dependent absorption;
- heterogeneous gray matter, white matter, nuclei, glia, and extracellular space;
- anisotropic or direction-dependent scattering;
- refractive-index boundaries at ACSF/slice/glass;
- lateral blur and photons entering from outside the nominal circular aperture;
- cell-specific scattering by the reconstructed PT5B neuron.

The single-cell morphology is used to sample an imposed irradiance field; it is not used to calculate tissue scattering.

## 6. ChR2 Conductance Model

### Current equation

ChR2 is a conductance, not a prescribed current source:

```text
I_ChR2 = gbar * open_fraction(t) * (V - E_ChR2)
E_ChR2 = 0 mV
```

`E_ChR2 = 0 mV` is the standard approximation for a nonselective cation conductance. It means inward current is strong at -70 mV and decreases as the membrane approaches 0 mV.

References:

- Nikolic K, et al. Photocycles of channelrhodopsin-2. Photochemistry and Photobiology. 2009. https://doi.org/10.1111/j.1751-1097.2008.00460.x
- Foutz TJ, Arlow RL, McIntyre CC. Theoretical principles underlying optical stimulation of a channelrhodopsin-2 positive pyramidal neuron. Journal of Neurophysiology. 2012. https://doi.org/10.1152/jn.00501.2011

### Four-state photocycle

The retained model has four principal states:

```text
C1: dark-adapted closed state
O1: primary open state
O2: light-adapted/desensitized open state
C2: light-adapted closed state
```

The open conductance is:

```text
open = O1 + gamma * O2
```

This captures light-driven opening, desensitization, a lower-conductance adapted open state, and dark recovery. It does not represent every spectroscopic intermediate of the ChR2 photocycle.

### Initial state

The default initial state is fully dark adapted:

```text
C1 = 1
O1 = O2 = C2 = 0
```

Wang does not report enough information about pre-illumination, YFP observation, pulse ordering, dark-adaptation intervals, or inter-trial intervals to reconstruct protocol history. A `C2_initial` sensitivity analysis was performed but was not adopted as the default.

### Photon-flux conversion

The mechanism converts irradiance into photon flux using:

```text
E_photon = h*c/lambda
photon_flux = irradiance / E_photon
```

and then maps photon flux to light-driven transition rates using retinal cross-section, loss, quantum-efficiency-like factors, and a finite activation variable.

Current implementation parameters include:

```text
wavelength = 470 nm
sigma_retinal = 12e-20 m2
wloss = 1.3
tau_chr2 = 1.3 ms
epsilon1 = 0.8535
epsilon2 = 0.14
```

These are transferred model parameters, not measurements from Thy1-ChR2 line 18. Their provenance is weaker than the Wang current targets and they remain a model-form uncertainty.

### Temperature scaling

All NEURON mechanisms are run with:

```text
h.celsius = 22 C
```

For ChR2, the final four-state branch uses parameter-specific Q10 factors previously extracted from a Williams H134R implementation:

```text
Q10 epsilon1 = 1.46
Q10 epsilon2 = 2.77
Q10 gd1 = 1.97
Q10 gd2 = 1.77
Q10 e12 = 1.10
Q10 e21 = 1.95
Q10 gr = 2.56
```

with:

```text
rate_22 = rate_37 / Q10^((37 - 22)/10)
```

This is not a direct wild-type ChR2 line-18 temperature measurement. It transfers H134R-family temperature dependence onto a wild-type Foutz/Nikolic-style scaffold. It is one of the largest remaining uncertainties, especially for sustained current and recovery.

The local repository does not currently contain a complete primary-paper provenance chain for every one of these numerical Q10 values. They should therefore be treated as implementation-derived priors, not as independently verified line-18 constants.

### Membrane localization

The final reference assumes uniform `gbar` density over illuminated retained soma and dendritic membrane. Axon is excluded.

Soma-only, dendrite-only, proximal-enriched, and proximal-reduced cases were tested as sensitivity analyses. No line-18 subcellular quantitative expression data were available to select one of them, so uniform density was retained as the least-committal assumption.

## 7. Meaning and Fitting of gbar

`gbar` is the maximum ChR2 conductance density in `mS/cm2`.

It absorbs several inseparable biological quantities:

- membrane ChR2 molecule density;
- functional fraction of expressed ChR2;
- single-channel conductance;
- membrane trafficking/localization;
- residual light-delivery calibration error;
- mismatch between the PT5B proxy and Wang's recorded cell.

### Calibration target

The final branch fits `gbar` to Wang's narrower saturating-current estimate:

```text
Imax = 0.642 nA
```

The numerical calibration pulse is:

```text
100 mW/mm2 for 100 ms
```

This `100 mW/mm2` pulse is a project-specific numerical saturation convention. Wang did not report this exact calibration pulse. It is used only to approach the asymptote; it is excluded from the seven-point intensity K fit.

After fitting, the same `gbar` is tested against Wang's 9.2 mW/mm2, 1 s trace and the seven-point intensity response.

For the soma-centered field the resulting effective value is approximately:

```text
gbar = 0.1009 mS/cm2
```

This is not a direct molecule count or universal Thy1-ChR2 line-18 constant.

## 8. Voltage Clamp and Recording Chain

### Somatic voltage clamp

The biological membrane model is clamped using a somatic NEURON `SEClamp`:

```text
command potential = -70 mV
series resistance = 10 MOhm
```

Wang reports pipette resistance `2-7 MOhm`, not access resistance after break-in. The adopted `10 MOhm` is a plausible effective access/series-resistance assumption, selected after sensitivity tests. It is not a measured value from Wang.

Because ChR2 current is distributed over a large dendritic tree, a somatic clamp does not perfectly hold all dendrites at -70 mV. This space-clamp limitation is retained rather than forcing `Rs = 0`.

General references for patch clamp and distributed voltage-clamp limitations:

- Sakmann B, Neher E. Patch clamp techniques for studying ionic channels in excitable membranes. Annual Review of Physiology. 1984. https://doi.org/10.1146/annurev.ph.46.030184.002323
- Li S, et al. Determination of effective synaptic conductances using somatic voltage clamp. PLOS Computational Biology. 2019. https://doi.org/10.1371/journal.pcbi.1006871

### Electrode/amplifier readout filter

After the raw `SEClamp` current is calculated, a phenomenological readout filter is applied:

```text
effective electrode resistance = 10 MOhm
effective capacitance = 100 pF
effective RC tau = R*C = 1.0 ms
then 4 cascaded first-order low-pass stages
tau per stage = 0.25 ms
```

This suppresses the very fast onset component and delays the apparent current peak.

Critical interpretation:

- `100 pF` should be viewed as a lumped effective recording-chain capacitance, not a measured literal glass-pipette capacitance from Wang;
- the four-stage `0.25 ms` filter is not an Axoclamp 2B setting reported by Wang;
- the values were chosen as a plausible phenomenological approximation and validated by their effect on the observable trace;
- this is not a full pipette-electrode feedback circuit or an exact amplifier model.

The filter therefore has weaker evidential status than the experimental protocol and should not be used to claim molecular ChR2 kinetics. General patch-clamp literature supports the existence of access-resistance, capacitance, bandwidth, and space-clamp distortions, but no cited paper uniquely supports the project's exact `100 pF`, `1 ms`, and four-stage filter values.

### Baseline subtraction

With original Dura-Bernal active biophysics, a matched no-light trace is subtracted from the light trace. This isolates the light-dependent clamp-current component and prevents intrinsic holding current from being interpreted as photocurrent.

## 9. Ionic Solutions and Concentrations

Wang's pipette and ACSF recipes were transcribed and used in a separate reversal-potential sensitivity analysis.

The final voltage-clamp calibration did not replace the full Dura-Bernal ion-concentration machinery with a dynamic Wang-solution model. It retains the Dura-Bernal membrane mechanisms and effective ionic settings.

No explicit modeling is included for:

- dynamic extracellular Na, K, Ca, or H accumulation/depletion;
- restricted extracellular space;
- pH shifts caused by ChR2 proton permeability;
- pipette dialysis over recording time;
- liquid-junction-potential dynamics.

This is acceptable for estimating a sub-nA photocurrent over short protocols, but it limits interpretation of repeated, prolonged, or very strong stimulation.

## 10. Synaptic and Network Inputs

The in vitro voltage-clamp model is a single isolated cell:

```text
no synaptic input
no recurrent network
no thalamic input
no neuromodulation
no spontaneous in vivo baseline state
```

This follows the intent of Wang's synaptic-blocker photocurrent measurements. It is not intended to reproduce in vivo excitability or motor output.

## 11. Numerical Assumptions

- NEURON multicompartment simulation.
- Original Dura-Bernal active mechanisms retained.
- Typical final integration time step: `0.1 ms` for field sweeps; some earlier reference runs used `0.05 ms`.
- ChR2 conductance installed only on selected illuminated segments.
- Irradiance is updated as an ideal rectangular pulse; the electronic shutter's measured temporal waveform is unavailable.
- The arc-lamp spectrum and shutter transient are not explicitly simulated.
- Hill K values are fitted from the Wang-like seven intensity points only.

## 12. Tested but Not Retained

The following were explored but are not part of the final reference assumption set:

- arbitrary `1x`/`4x` expression scales;
- the erroneous 500 um radial morphology cutoff;
- in vivo `mu_eff = 2.12 mm^-1` for the final slice model;
- fitting irradiance scale solely to force agreement;
- changing the ideal light rise time solely to fit time-to-peak;
- non-dark-adapted `C2_initial` as a default;
- ad hoc slow adaptation/desensitization branches;
- manually fitting `e12/e21/gd2/gamma` without independent evidence;
- PyRhO six-state results obtained before the morphology truncation issue was corrected;
- using Wang Fig. 3 firing as the primary `gbar` target;
- tuning PT5B leak/HCN solely to force Wang Table 1 input resistance and resting voltage.

## 13. What Is Well Constrained

Relatively strong constraints:

- Wang optical protocols and reported photocurrent summary values;
- room-temperature range;
- large-field area order;
- holding potential;
- line-18 expression and cortical L5 pyramidal-cell context;
- Dura-Bernal PT5B morphology and active membrane implementation;
- conductance-based rather than fixed-current ChR2 formulation.

## 14. Main Remaining Uncertainties

The largest unresolved uncertainties are:

1. Wang cell identity versus M1 PT5B proxy.
2. Exact arc-lamp spectrum, field shape, and tissue-plane irradiance calibration.
3. Effective in vitro optical attenuation and lateral scattering.
4. Wild-type ChR2 temperature scaling at 22 C.
5. Subcellular ChR2 localization in line 18.
6. Exact access resistance, capacitance compensation, amplifier bandwidth, and filtering used for Fig. 2.
7. Protocol history and dark adaptation.
8. Full electrical removal and sealing of dendrites cut by the slice.

## 15. Recommended Interpretation

The final `gbar ~= 0.10 mS/cm2` should be used as an order-of-magnitude, geometry-conditional effective ChR2 conductance density for the next in vivo modeling phase.

It is supported by agreement with Wang's peak-current scale and intensity K under a plausible soma-centered field. It is not validated as a unique molecular expression density because optical attenuation, cell identity, subcellular localization, recording filtering, and temperature scaling remain partially confounded.

The remaining mismatch in prolonged kinetics should be preserved as model uncertainty rather than eliminated through unsupported parameter fitting:

```text
model decay tau approximately 27 ms versus Wang approximately 48 ms
model duration K approximately 1.68 ms versus Wang approximately 3.2 ms
```

## Primary References

1. Wang H, et al. PNAS. 2007. https://doi.org/10.1073/pnas.0700384104
2. Dura-Bernal S, et al. Cell Reports. 2023. https://doi.org/10.1016/j.celrep.2023.112574
3. Yona G, et al. eNeuro. 2016. https://doi.org/10.1523/ENEURO.0059-15.2015
4. Nikolic K, et al. Photochemistry and Photobiology. 2009. https://doi.org/10.1111/j.1751-1097.2008.00460.x
5. Foutz TJ, et al. Journal of Neurophysiology. 2012. https://doi.org/10.1152/jn.00501.2011
6. Sakmann B, Neher E. Annual Review of Physiology. 1984. https://doi.org/10.1146/annurev.ph.46.030184.002323
7. Li S, et al. PLOS Computational Biology. 2019. https://doi.org/10.1371/journal.pcbi.1006871
