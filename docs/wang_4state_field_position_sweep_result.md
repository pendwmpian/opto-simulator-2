# Wang 4-state field-position sweep result

## Purpose

Test how strongly Wang et al. 2007-like voltage-clamp readouts depend on the assumed position of the `0.4 mm2` circular illumination field on the slice face.

This is a geometry sensitivity run, not a new ChR2 kinetics fit.

## Fixed Assumptions

- Cell: `PT5B_full`.
- Morphology: no radial morphology truncation. Axon excluded as before.
- Slice: Wang-style 300 um acute slice abstraction.
- Illumination: circular `0.4 mm2` field on the slice face.
- Field radius: `356.8 um`.
- Optical attenuation: in vitro `mu_eff = 1.3 mm^-1`.
- ChR2 model: four-state Foutz/Nikolic-style mechanism with Williams-style 37 -> 22 C Q10 scaling.
- Temperature: `22 C`.
- Voltage clamp: soma `SEClamp`, hold `-70 mV`.
- Electrode/readout: `Rs = 10 MOhm`, electrode Rs `10 MOhm`, pipette capacitance `100 pF`, amplifier filter tau `0.25 ms`.
- `gbar` calibration: refit at each field position to Wang's saturated photocurrent target, `Imax = 0.642 nA`, using `100 mW/mm2`, `100 ms`.

## Sweep

The default field center was the area-weighted apical dendrite centroid on the illumination plane.

Five field positions were tested:

- center
- `-200 um` along field axis-a
- `+200 um` along field axis-a
- `-200 um` along field axis-b
- `+200 um` along field axis-b

With the default `illumination_axis = z`, the plotted illumination plane is `x-y`.

## Outputs

- Summary CSV: `outputs/wang_4state_field_position_sweep/field_position_summary.csv`
- Intensity CSV: `outputs/wang_4state_field_position_sweep/field_position_intensity.csv`
- Duration CSV: `outputs/wang_4state_field_position_sweep/field_position_duration.csv`
- Summary figure: `outputs/wang_4state_field_position_sweep/field_position_sweep_summary.png`
- Per-position composite figures: `outputs/wang_4state_field_position_sweep/*_composite.png`

## Key Results

| position | gbar (mS/cm2) | 9.2 mW peak (nA) | steady/peak | TTP (ms) | tau (ms) | intensity K (mW/mm2) | duration K (ms) | illuminated area fraction | apic illuminated fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| apical center | 0.466 | 0.616 | 0.232 | 10.3 | 76.1 | 0.743 | 1.60 | 0.590 | 0.970 |
| axis-a -200 | 0.667 | 0.409 | 0.201 | 11.5 | 90.6 | 0.933 | 1.70 | 0.392 | 0.671 |
| axis-a +200 | 0.610 | 0.408 | 0.213 | 11.0 | 87.9 | 0.765 | 1.61 | 0.440 | 0.746 |
| axis-b -200 | 0.101 | 0.579 | 0.228 | 9.9 | 27.5 | 0.858 | 1.68 | 0.768 | 0.602 |
| axis-b +200 | 0.956 | 0.050 | 0.091 | 11.6 | 88.2 | 0.710 | 1.58 | 0.346 | 0.593 |

## Interpretation

The field position is not a harmless detail.

Moving the same `0.4 mm2` circle changes:

- fitted `gbar` by almost an order of magnitude,
- predicted 9.2 mW/mm2 peak current by more than 10-fold after saturated-current calibration,
- apparent time-to-peak by about `1.7 ms`,
- apparent inactivation tau by more than 3-fold,
- intensity K by about `0.71-0.93 mW/mm2`.

The axis-b `+200 um` position is a clear bad geometry for this morphology because it illuminates relatively little effective responsive membrane under the current slice/field convention. It can still be forced to match saturated current by increasing `gbar`, but then the 9.2 mW/mm2 response becomes far too small.

The most useful comparison set is probably:

- apical center,
- axis-a `-200`,
- axis-a `+200`,
- axis-b `-200`.

The apical-center and axis-b `-200` positions both produce Wang-like intensity K values, but they imply very different `gbar` and decay tau. This means field geometry should be constrained before using Wang Fig. 2 timing or K values to tune ChR2 kinetics.

## Caveat

This run uses the restored full morphology with no radius truncation. Earlier PyRhO6 sensitivity results that used default `truncate_radius_um = 500` should not be interpreted as final geometry-controlled Wang calibration results.
