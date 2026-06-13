# Full-model quiet PT5B baseline: run report

Updated: 2026-06-13 UTC

## Objective

Run the Dura-Bernal et al. (2023) full NetPyNE/NEURON network while retaining
natural network dynamics, all-cell spikes, and one PT5B soma Vm trace. Avoid
OOM during connection creation by reducing recorded and serialized metadata,
not by replacing or pruning the network model.

## Model and environment

- Upstream repository: `external/M1_NetPyNE_CellReports_2023`
- Upstream commit: `f8ed26418a3120b2cb48dacebf8a96afb6afa9a9`
- Network scale: upstream `0.3` configuration (12,173 cells)
- MPI ranks: 8, bound to hardware threads
- Python: 3.13.14
- NetPyNE: 1.1.1
- NEURON: 9.0.1, source-built with OpenMPI support
- Trial 0 seeds: conn 4321, stim 1234, loc 4321
- Resource limits: `MemoryHigh=200G`, `MemoryMax=220G`
- Thread limits: OMP/OpenBLAS/MKL each set to 1 thread

## Run-only configuration

Retained:

- Full cell models, morphology, ion channels, synaptic mechanisms and NetCons
- Background network inputs and stimulation mechanisms
- ChR2 mechanisms in the compiled mechanism set
- Spikes from all cells
- Soma Vm from the lowest-GID PT5B cell, GID 5130

Disabled or omitted from saved output:

- Vm and state traces from all non-target cells
- LFP, dipole and stimulus traces
- Cell section and connection metadata
- Serialized network and NetParams
- In-run analysis and plotting

The target cell is owned by MPI rank 2. Only rank 2 records one Vm trace; all
other ranks report zero traces. Internal states and synapses remain active in
the simulator.

## Completed validation runs

### Full-network build-only

- Directory: `full-model-runs/full-build-trial0`
- Status: success
- Wall time: 12 min 21 sec (07:36:17--07:48:38 UTC)
- Cell creation: about 31 sec
- Connection creation and subcellular placement: about 11 min 21 sec
- Peak monitored aggregate RSS: 78.60 GB
- Peak single-rank RSS after connection creation: about 9.90 GB
- Swap: none

This confirms that full connection creation completes below the configured
memory limits. The earlier pueue task 0 failed before model startup because
OpenMPI counted only four physical-core slots. The successful command uses
`--use-hwthread-cpus --bind-to hwthread`; this is an execution correction, not
a model change.

### Full-network 10 ms dynamics

- Directory: `full-model-runs/full-dynamics-10ms-trial0`
- Status: success
- Wall time: 16 min 30 sec (07:49:32--08:06:02 UTC)
- Approximate integration time after construction: 3 min 30 sec
- Peak monitored aggregate RSS: 81.08 GB
- Target Vm samples: 100 at 0.1 ms spacing, 0.0--9.9 ms
- Target Vm range: -77.046 to -70.043 mV
- Non-finite target Vm values: none
- All-network spikes in the first 10 ms: 0

Zero spikes over this short initialization interval are not sufficient for a
baseline firing-rate conclusion. The completed run verifies MPI dynamics,
minimal recording, gathering, and output serialization.

## Completed 1 s baseline run

- Directory: `full-model-runs/full-quiet-1s-trial0`
- pueue simulation task: 5 (`full-quiet-1s-trial0`), success
- pueue monitor task: 6 (`monitor-full-quiet-1s-trial0`), success
- Start: 2026-06-12 08:10:46 UTC
- End: 2026-06-12 14:09:34 UTC
- Requested biological duration: 1000 ms
- Wall time: 5 h 58 min 47 sec
- All-network spikes: 46,802 from 5,869 / 12,173 cells
- Mean rate over all modeled cells: 3.845 Hz
- Target PT5B GID 5130 spikes: 0
- Target Vm mean after the first 100 ms: -66.975 mV
- Target Vm SD after the first 100 ms: 2.272 mV
- Target Vm range over the full run: -76.478 to -61.472 mV
- Peak monitored aggregate RSS: 81.47 GB
- Swap: none

## Corrected scientific interpretation

This run is a successful full-biophysics, 0.3-scale network execution and a
successful test of minimal recording and OOM avoidance. It is **not** a
reproduction of Schiemann et al. (2015) quiet wakefulness.

Schiemann et al. reported an L5B mean Vm of -51.1 +/- 0.8 mV during quiet
wakefulness. Identified PT-type cells had mean Vm -47.6 +/- 7.1 mV, Vm SD
3.3 +/- 0.8 mV, and firing rate 6.7 +/- 4.7 Hz (n=6). The selected model cell
was therefore about 19 mV too hyperpolarized and did not fire.

The executed condition also differed from the upstream Figure 2 representative
quiet condition: it retained `ihGbar=1.0` instead of 0.75, simulated only 1 s
instead of 5 s, included initialization transient, and used one seed. Upstream
Figure 2 evaluates 1--5 s over 5 connectivity x 5 stimulation seeds. In
addition, GID 5130 was selected only because it was the lowest PT5B GID, not
because it was physiologically representative.

The correct conclusion is:

- MPI full-network execution: successful
- OOM avoidance and minimal output: successful
- Dura-Bernal Figure 2 quiet reproduction: not yet tested correctly
- Schiemann quiet-wake L5B/PT physiology: not reproduced in this trial

## Assumptions and interpretation limits

- This is the normal/control quiet condition; explicit pulse, IClamp and
  NetStim test stimulation are disabled.
- ChR2 is compiled and available but is not illuminated in this baseline run.
- `scale=0.3` is retained from the upstream runnable configuration and is not a
  one-cell or surrogate reduction.
- A corrected run should use `ihGbar=0.75`, at least 5 s duration, a 1--5 s
  analysis window, and PT5B population statistics. Multiple Vm targets should
  be inspected before fixing one representative cell for the optical stage.
- Trial 1 and trial 2 should only be launched after trial 0 output and runtime
  are reviewed, to avoid spending two additional multi-hour runs on a faulty
  condition.

## Reproduction command

```bash
mpiexec --use-hwthread-cpus --bind-to hwthread --mca btl self,vader,tcp -n 8 \
  .venv/bin/python scripts/run_full_quiet_pt5b.py \
  --trial 0 --duration-ms 1000 \
  --output-dir docs-temp/full-model-runs/full-quiet-1s-trial0
```
