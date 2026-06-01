# Wang solution-derived reversal potential check

## Goal

Re-check Wang et al. 2007 slice patch solution conditions, compute approximate Nernst reversal potentials from those solutions, and test whether changing only ionic reversal potentials can explain the discrepancy between the Dura-Bernal PT5B/IT5B intrinsic responses and Wang Fig. 3/Table 1.

## Wang et al. 2007 solution conditions

Source: Wang et al. 2007, PNAS, Methods.

Recording pipette:

- 130 mM K-gluconate
- 2 mM NaCl
- 4 mM MgCl2
- 20 mM HEPES
- 4 mM Na2ATP
- 0.4 mM NaGTP
- 250 uM K-EGTA
- pH 7.3

Extracellular solution:

- 125 mM NaCl
- 2.5 mM KCl
- 2 mM CaCl2
- 1 or 1.3 mM MgCl2
- 20 mM dextrose/glucose
- 1.25 mM NaH2PO4
- 26 mM NaHCO3
- pH 7.4 after bubbling with 95% O2 / 5% CO2

Wang reference physiology:

- ChR2+ resting Vm: -58.3 +/- 1.8 mV
- ChR2+ input resistance: 126 +/- 13 MOhm
- ChR2- resting Vm: -61.4 +/- 1.4 mV
- ChR2- input resistance: 110 +/- 4.5 MOhm

## Reversal-potential assumptions

Use the Nernst equation at 22 deg C:

`E = (RT / zF) ln([out] / [in])`

Important caveat: the pipette recipe is a total chemical recipe, not a direct measurement of free intracellular ion activity after dialysis. Therefore Na has two useful bounds:

- `wang_total_na`: internal Na includes NaCl + Na2ATP + NaGTP sodium, i.e. 2 + 8 + 0.8 = 10.8 mM.
- `wang_free_na2`: internal Na is treated as the 2 mM NaCl term only. This is an upper-bound E_Na sensitivity check, not necessarily the most physiological value.

For K, use 130 mM K-gluconate plus 0.25 mM K-EGTA as a small total K correction.

For Ca, the free intracellular Ca concentration is not specified by the recipe because EGTA buffers it. Use 100 nM as a standard free-Ca nominal value for the main Wang-derived profile, and also keep a `wang_total_na_dura_ca` variant where E_Ca is left near the Dura-Bernal effective value.

## Simulation

Keep these fixed:

- Dura-Bernal morphology and channel densities.
- Original Dura-Bernal mechanisms.
- No ChR2.
- 22 deg C.
- The existing slice pruning abstraction with soma entry depth 25 um, because that was the strongest PT5B Rin case.

Sweep these ion profiles:

- `dura_default`: no reversal-potential override.
- `wang_total_na`: E_Na from total pipette Na, E_K from Wang K, E_Ca from 100 nM free Ca.
- `wang_free_na2`: E_Na from 2 mM internal Na, E_K from Wang K, E_Ca from 100 nM free Ca.
- `wang_total_na_dura_ca`: E_Na and E_K from Wang solution, E_Ca left at Dura effective value.

Readouts:

- Resting Vm.
- Apparent input resistance from a 50 pA step.
- Spike count and latency for 1-s somatic current steps.
- AP overshoot/max Vm as a sanity check.

