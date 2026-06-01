# Wang solution-derived reversal potential sweep result

## Source conditions

Wang et al. 2007 reports the following whole-cell patch conditions for line 18 slice recordings:

- Pipette: 130 mM K-gluconate, 2 mM NaCl, 4 mM MgCl2, 20 mM HEPES, 4 mM Na2ATP, 0.4 mM NaGTP, 250 uM K-EGTA, pH 7.3.
- ACSF: 125 mM NaCl, 2.5 mM KCl, 2 mM CaCl2, 1 or 1.3 mM MgCl2, 20 mM dextrose/glucose, 1.25 mM NaH2PO4, 26 mM NaHCO3, pH 7.4 with 95% O2 / 5% CO2.
- Pipette resistance: 2-7 MOhm.
- Junction potential: 10 mV taken into account for reported membrane potentials.
- Temperature: 21-24 deg C.

## Computed Nernst values at 22 deg C

The pipette recipe does not directly specify free intracellular Na or free intracellular Ca activity, so two Na interpretations were tested.

| profile | E_Na (mV) | E_K (mV) | E_Ca (mV) | assumption |
|---|---:|---:|---:|---|
| dura_default | Dura | Dura | Dura | original Dura-Bernal ion settings |
| wang_total_na | 67.30 | -100.55 | 125.94 | Na_in = 2 + 8 + 0.8 mM; Ca_in = 100 nM |
| wang_free_na2 | 110.19 | -100.55 | 125.94 | Na_in = 2 mM only; high E_Na sensitivity check |
| wang_total_na_dura_ca | 67.30 | -100.55 | Dura | Na/K from Wang, Ca left unchanged |

Notes:

- Na_out was treated as 125 + 26 + 1.25 = 152.25 mM.
- K_in was treated as 130 + 0.25 = 130.25 mM.
- Ca_in is not directly known from the recipe because EGTA buffers free Ca. 100 nM is a nominal free-Ca choice, not a measured Wang value.

## Simulation setup

- Cell models: Dura-Bernal PT5B full and IT5B full.
- Biophysics: original Dura-Bernal mechanisms and channel densities.
- ChR2: absent.
- Temperature: 22 deg C.
- Slice pruning: enabled, soma entry depth 25 um, slice axis z, truncate radius 500 um.
- Perturbation: only section-level `ena`, `ek`, and `eca` were overridden where that ion existed.

Outputs:

- `outputs/wang_solution_reversal_sweep_pt5b/`
- `outputs/wang_solution_reversal_sweep_it5b/`

## PT5B result

| profile | rest Vm (mV) | Rin (MOhm) | 0.2 nA | 0.3 nA | 0.5 nA | main effect |
|---|---:|---:|---:|---:|---:|---|
| dura_default | -75.70 | 72.52 | 0 Hz | 15 Hz | 30 Hz | baseline |
| wang_total_na | -75.66 | 72.86 | 0 Hz | 14 Hz | 24 Hz | larger AP overshoot, slightly lower f-I |
| wang_free_na2 | -75.66 | 72.92 | 1 Hz | 14 Hz | 22 Hz | very large AP overshoot, modest f-I change |
| wang_total_na_dura_ca | -75.66 | 72.86 | 0 Hz | 14 Hz | 24 Hz | same as wang_total_na |

PT5B conclusion:

- Matching Wang solution-derived reversal potentials alone does not move PT5B toward Wang Table 1 resting Vm or Rin.
- Resting Vm remains near -75.7 mV, far from Wang ChR2+ -58.3 mV.
- Rin remains near 73 MOhm, far below Wang ChR2+ 126 MOhm and ChR2- 110 MOhm.
- Increasing E_Na mainly increases AP overshoot and modestly changes the f-I curve. It does not explain the resting-state mismatch.

## IT5B result

| profile | rest Vm (mV) | Rin (MOhm) | 0.1 nA | 0.2 nA | 0.3 nA | main effect |
|---|---:|---:|---:|---:|---:|---|
| dura_default | -80.39 | 237.11 | 8 Hz | 17 Hz | 24 Hz | baseline |
| wang_total_na | -80.13 | 243.76 | 8 Hz | 16 Hz | 22 Hz | larger AP overshoot, slightly lower f-I |
| wang_free_na2 | -80.13 | 246.17 | 9 Hz | 16 Hz | 22 Hz | very large AP overshoot |
| wang_total_na_dura_ca | -80.13 | 243.76 | 8 Hz | 16 Hz | 22 Hz | same as wang_total_na |

IT5B conclusion:

- IT5B remains much more excitable than PT5B and gives Wang-like firing frequency around 0.3 nA.
- However, its Rin is too high and rest Vm is too hyperpolarized relative to Wang Table 1.
- As with PT5B, solution-derived reversals do not solve the resting-state mismatch.

## Interpretation

This result argues against the idea that the Wang-vs-Dura discrepancy is mainly caused by using the wrong external/internal solution reversal potentials. The solution correction is still important for AP shape and ionic driving force, but the large mismatch in rest Vm and Rin likely requires differences in one or more of:

- cell class or developmental state of Wang's recorded cortical L5 neurons,
- passive leak/resting conductance,
- HCN and other subthreshold active conductances,
- morphology retained in the acute slice,
- recording selection/state not captured by the Dura-Bernal in vivo network cell parameters.

For ChR2 expression-order calibration, this means Wang Fig. 2 photocurrent can still constrain a line18 membrane conductance scale, but Wang Fig. 3 current-clamp excitability should not be treated as a direct PT5B/M1 constraint unless the intrinsic cell state is separately matched.

