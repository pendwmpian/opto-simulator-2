# PT5B Slice Position Input-Resistance Sweep

## Purpose

Estimate where the Dura-Bernal et al. 2023 M1 `PT5B_full` morphology should sit inside a Wang et al. 2007-like acute slice before adding ChR2.

The immediate target is the ChR2-negative / no-ChR2 electrical baseline. The slice placement is judged by somatic input resistance, with Wang et al. 2007's reported target set to 126 MOhm.

## Model Scope

- Cell: `external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl`
- Biophysics: default simplified biophysics from `CellFromNetPyNE`, unless `--use-original-biophysics` is passed.
- ChR2: not installed.
- Slice: 300 um Wang-style vertical acute slice, using `segment_rows_vertical_slice`.
- Slice axes: `z` and `x` are swept as candidate slice-thickness axes.
- Off-slice handling: density-mechanism conductances are set to zero for off-slice segments; sections with no retained segments get near-zero capacitance.
- Rin measurement: somatic hyperpolarizing current step, default `-0.05 nA`; steady-state voltage deflection divided by current amplitude. Since mV/nA equals MOhm, no additional unit conversion is needed.
- Optional leak scale: `--g-pas-scale` scales both `pas.g` and the simplified model's `hh.gl`, because both contribute to baseline membrane leak.

## Experiment Command

```bash
rtk .venv/bin/python scripts/sweep_pt5b_slice_rin.py
```

The default sweep tests soma entry depths:

```text
25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275 um
```

for both `z` and `x` slice-thickness axes.

## Outputs

- `outputs/pt5b_slice_rin_sweep/pt5b_slice_rin_sweep.csv`
- `outputs/pt5b_slice_rin_sweep/best_pt5b_slice_rin.json`
- `outputs/pt5b_slice_rin_sweep/best_trace.csv`
- `outputs/pt5b_slice_rin_sweep/pt5b_slice_rin_sweep.png`

## Initial Result

Default simplified leak (`--g-pas-scale 1.0`) could not reach the Wang target by slice position alone.

Best default geometry-only result:

- output directory: `outputs/pt5b_slice_rin_sweep`
- slice thickness axis: `x`
- soma entry depth: `25 um`
- input resistance: `70.69 MOhm`
- target error: `55.31 MOhm`
- retained membrane area fraction: `0.464`

This indicates that, with the current simplified PT5B membrane leak, slice placement alone is insufficient to reproduce the 126 MOhm ChR2-negative baseline.

## Leak-Adjusted Follow-Up

Because the simplified model includes both `pas.g` and `hh.gl_hh` leak terms, a follow-up sweep scaled both by `--g-pas-scale`.

Coarse follow-up:

```bash
rtk .venv/bin/python scripts/sweep_pt5b_slice_rin.py --out-dir outputs/pt5b_slice_rin_sweep_leak0p5 --g-pas-scale 0.5
```

Best coarse result:

- slice thickness axis: `z`
- soma entry depth: `50 um`
- input resistance: `123.93 MOhm`
- target error: `2.07 MOhm`

Refined follow-up:

```bash
rtk .venv/bin/python scripts/sweep_pt5b_slice_rin.py --out-dir outputs/pt5b_slice_rin_sweep_leak0p49_refined --g-pas-scale 0.49 --slice-axes z,x --soma-entry-depths-um 25,40,50,60,75,250,275
rtk .venv/bin/python scripts/sweep_pt5b_slice_rin.py --out-dir outputs/pt5b_slice_rin_sweep_leak0p485_refined --g-pas-scale 0.485 --slice-axes z,x --soma-entry-depths-um 25,40,50,60,75,250,275
```

Best refined result:

- output directory: `outputs/pt5b_slice_rin_sweep_leak0p485_refined`
- slice thickness axis: `x`
- soma entry depth: `60 um`
- slice center: `90.15 um`
- current step: `-0.05 nA`
- baseline voltage: `-73.03 mV`
- steady voltage: `-79.32 mV`
- input resistance: `125.90 MOhm`
- target error: `0.10 MOhm`
- retained membrane area: `18226.80 um2`
- total membrane area: `28919.36 um2`
- retained membrane area fraction: `0.630`

## Interpretation Plan

1. Treat `x` slice thickness axis, `60 um` soma entry depth, and `g_pas_scale = 0.485` as the current ChR2-negative baseline candidate.
2. Keep the geometry-only result visible: it shows that this morphology/slice placement does not explain the input-resistance mismatch without also reducing leak.
3. Use the selected no-ChR2 baseline position for the next Wang current-clamp / voltage-clamp ChR2-positive simulations.
