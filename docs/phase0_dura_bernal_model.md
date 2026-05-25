# Phase 0: Dura-Bernal 2023 M1 Model Intake

## Source

Downloaded GitHub repository:

```text
https://github.com/suny-downstate-medical-center/M1_NetPyNE_CellReports_2023
```

The repository README describes this as a multiscale mouse primary motor cortex model developed with NetPyNE and NEURON. It also states that the full model requires approximately 2 hours on 96 cores for 1 second of simulation.

## Local Layout

The repository is stored locally under:

```text
external/M1_NetPyNE_CellReports_2023
```

Important directories:

```text
sim/netParams.py   network and cell-rule assembly
sim/cfg.py         simulation configuration
sim/cells/         cell rules and cell source files
sim/conn/          connectivity and weight normalization files
sim/mod/           NMODL ion-channel and synapse mechanisms
data/conn/         connectivity and source cellParams used by analysis
```

## What Form Are The Cell Models In?

The model contains both full and reduced L5B-related cell rules.

The most directly useful L5B candidates are:

```text
sim/cells/PT5B_full_cellParams.pkl
sim/cells/PT5B_reduced_cellParams.pkl
sim/cells/IT5B_full_cellParams.pkl
sim/cells/IT5B_reduced_cellParams.pkl
```

These `.pkl` files are NetPyNE cell parameter rules. They are not just scalar parameter tables. For the full cells, each section includes 3D morphology points in `geom["pt3d"]`, so the cell shape can be reconstructed and visualized without rerunning the full network.

The repository also contains source mechanisms:

```text
sim/cells/PTcell.hoc
sim/cells/PTcell.py
sim/cells/ITcell.py
sim/cells/CSTR6.py
sim/cells/SPI6.py
sim/mod/*.mod
```

For simulation, the NMODL mechanisms in `sim/mod` must be compiled before running detailed NEURON dynamics.

## L5B Candidate Summary

### PT5B Full

```text
file: sim/cells/PT5B_full_cellParams.pkl
sections: 174
section counts: apic 103, dend 69, soma 1, axon 1
sections with pt3d: 174
pt3d points: 2620
total section length: 11552.72 um
mechanisms: cadad, cal, can, hd, kBK, kap, kdmc, kdr, nax, pas, savedist
```

This is the best first candidate for a large L5B pyramidal tract style cell with rich apical dendritic morphology.

### IT5B Full

```text
file: sim/cells/IT5B_full_cellParams.pkl
sections: 95
section counts: apic 27, dend 54, soma 1, axon 13
sections with pt3d: 95
pt3d points: 2723
total section length: 6859.54 um
mechanisms: cadad, cal, can, cat, ican, ih, kBK, kap, kdr, nap, nax, pas
```

This is the best first candidate for a L5B intratelencephalic pyramidal cell.

### Reduced Cells

Both `PT5B_reduced` and `IT5B_reduced` are 6-section reduced models:

```text
soma, Adend1, Adend2, Adend3, Bdend, axon
```

They do include simple `pt3d` coordinates, but these are schematic rather than reconstructed dendritic trees. They are useful for quick ChR2 mechanism development and parameter sweeps, but not ideal for spatial receptive-field estimation.

## Recommended First Extraction

Use two tracks:

1. Fast development track: `PT5B_reduced` or `IT5B_reduced`
2. Spatially meaningful track: `PT5B_full` first, then `IT5B_full`

For the optogenetic question, `PT5B_full` is the most valuable first detailed target because it has extensive apical and basal dendritic morphology and is a canonical L5B output-like pyramidal population. `IT5B_full` should be kept as a second candidate because M1 L5B responses in the source model distinguish IT5B and PT5B populations.

## Generated Visualizations

The morphology visualizations were generated from the stored NetPyNE `pt3d` geometry.

```text
outputs/figures/PT5B_full_morphology.png
outputs/figures/IT5B_full_morphology.png
```

Color coding:

```text
black: soma
green: axon
red: apical dendrite
blue: basal dendrite
```

## Current Phase 0 Status

Completed:

- Created Python 3.13 virtual environment with `uv venv --python=3.13`.
- Installed `netpyne`, `neuron`, `numpy`, and `matplotlib` into the virtual environment.
- Downloaded the Dura-Bernal 2023 M1 model repository.
- Confirmed that full L5B candidate cells contain reconstructable 3D morphology.
- Generated 3D/projection morphology figures for PT5B full and IT5B full.

Not yet done:

- Compile the model's NMODL mechanisms with `nrnivmodl`.
- Instantiate the cells inside NEURON and run electrophysiology validation.
- Add ChR2 mechanism.

