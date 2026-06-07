# Wang 4-state soma-centered field fine sweep result

## Purpose

Repeat the Wang Fig. 2-like field-position sweep under the more likely experimental assumption that the wide-field illumination was centered near the patched soma.

This is a narrower geometry sensitivity test than the previous apical-centered sweep.

## Fixed Assumptions

- Cell: `PT5B_full`.
- Morphology: no radial morphology truncation. Axon excluded.
- Slice: Wang-style 300 um acute slice abstraction.
- Illumination: circular `0.4 mm2` field on the slice face.
- Field radius: `356.8 um`.
- Field center: soma-centered, swept in a `3 x 3` grid.
- Offset grid: `-100, 0, +100 um` along each field-plane axis.
- Optical attenuation: in vitro `mu_eff = 1.3 mm^-1`.
- ChR2 model: four-state Foutz/Nikolic-style mechanism with Williams-style 37 -> 22 C Q10 scaling.
- Temperature: `22 C`.
- Voltage clamp/readout: same settings as the previous field-position sweep.
- `gbar` calibration: refit at each field position to Wang's saturated photocurrent target, `Imax = 0.642 nA`, using `100 mW/mm2`, `100 ms`.

## Outputs

- Summary CSV: `outputs/wang_4state_field_position_soma_fine_sweep/field_position_summary.csv`
- Intensity CSV: `outputs/wang_4state_field_position_soma_fine_sweep/field_position_intensity.csv`
- Duration CSV: `outputs/wang_4state_field_position_soma_fine_sweep/field_position_duration.csv`
- Summary figure: `outputs/wang_4state_field_position_soma_fine_sweep/field_position_sweep_summary.png`
- Per-position composite figures: `outputs/wang_4state_field_position_soma_fine_sweep/*_composite.png`

## Key Results

Across the soma-centered `3 x 3` sweep:

- `gbar`: essentially constant at `0.1009 mS/cm2`
- 9.2 mW/mm2 1 s peak current: `0.573-0.578 nA`
- steady/peak: `0.228-0.229`
- time-to-peak: `9.9 ms`
- inactivation tau: `27.0-27.5 ms`
- intensity K: `0.858-0.859 mW/mm2`
- duration K: `1.680 ms`
- illuminated membrane fraction: `0.64-0.74`
- illuminated apical fraction: `0.39-0.55`

## Interpretation

The soma-centered field is locally robust. Moving a `0.4 mm2` circular field by `+/-100 um` around the patched soma barely changes the fitted or predicted readouts.

This is qualitatively different from the broader apical-centered sweep, where moving the field by `200 um` could change effective illuminated membrane area and fitted `gbar` substantially.

For Wang-style patch-clamp calibration, the soma-centered geometry is therefore a more stable default assumption:

- It is experimentally plausible because the operator would normally illuminate around the patched soma unless stated otherwise.
- It avoids large arbitrary changes in `gbar` caused only by field placement.
- It gives an intensity K close to Wang's reported `0.84 mW/mm2`.

However, it still leaves two kinetic mismatches:

- decay tau remains too fast: about `27 ms` versus Wang's reported `48 ms`;
- duration K remains too low: about `1.68 ms` versus Wang's reported `3.2 ms`.

These remaining mismatches should be treated as ChR2 kinetics / recording-filter / protocol-history issues, not as field-position artifacts within a soma-centered `0.4 mm2` illumination model.
