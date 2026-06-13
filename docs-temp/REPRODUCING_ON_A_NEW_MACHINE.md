# Reproducing the full-network runs

## Repository state

The runner expects the Dura-Bernal model at:

`external/M1_NetPyNE_CellReports_2023`

Clone that repository and check out commit:

`f8ed26418a3120b2cb48dacebf8a96afb6afa9a9`

```bash
git clone https://github.com/suny-downstate-medical-center/M1_NetPyNE_CellReports_2023.git \
  external/M1_NetPyNE_CellReports_2023
git -C external/M1_NetPyNE_CellReports_2023 checkout \
  f8ed26418a3120b2cb48dacebf8a96afb6afa9a9
```

The `external/` directory and `.venv/` are intentionally not committed.

## Python and MPI

The tested environment used Python 3.13.14, OpenMPI 4.1.6, mpi4py 4.1.2,
NetPyNE 1.1.1, and NEURON 9.0.1. NEURON must be built with MPI enabled; the
stock wheel may not provide the required MPI configuration.

Create the uv environment from the committed project files, then rebuild
mpi4py and NEURON against the machine's OpenMPI installation. The tested
NEURON source configuration was:

```bash
CMAKE_ARGS='-DNRN_ENABLE_MPI=ON -DNRN_ENABLE_MPI_DYNAMIC=OFF -DNRN_ENABLE_INTERVIEWS=OFF' \
  uv pip install --python .venv/bin/python --no-binary neuron neuron==9.0.1
```

Compile the upstream MOD files from its `sim/` directory using the resulting
`nrnivmodl`. Confirm that a two-rank smoke test distributes cells before
launching the full model.

## Tested command

```bash
mpiexec --use-hwthread-cpus --bind-to hwthread --mca btl self,vader,tcp -n 8 \
  .venv/bin/python scripts/run_full_quiet_pt5b.py \
  --trial 0 --duration-ms 1000 \
  --output-dir docs-temp/full-model-runs/full-quiet-1s-trial0
```

Use pueue and a memory-controlled systemd scope for long runs, as recorded in
`full_model_quiet_pt5b_run_report.md`.

## Important scientific correction

The completed 1 s run is a technical baseline only. It is not the upstream
Figure 2 quiet condition and does not reproduce Schiemann et al. (2015). A
corrected comparison should set `ihGbar=0.75`, run at least 5 s, discard the
first second, and evaluate PT5B population statistics over 1--5 s.
