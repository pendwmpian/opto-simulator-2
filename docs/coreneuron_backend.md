# CoreNEURON backend

The full-network runner now defaults to the CPU CoreNEURON backend. Build the
model mechanisms first:

```bash
uv run python scripts/build_coreneuron_mechanisms.py
```

The build creates a temporary merged MOD source directory under the ignored
upstream checkout, applies thread-safety compatibility changes, installs a
CoreNEURON-compatible VecStim, includes the project ChR2 mechanisms, and runs
`nrnivmodl -coreneuron`.

CoreNEURON 8.2 cannot execute the per-timestep Python callbacks used by NetPyNE
1.1.1 for LFP calculation. The runner therefore rejects CoreNEURON plus LFP
instead of writing a plausible-looking but invalid array. Use
`--no-record-lfp`. Soma Vm is registered at the integration step and then
downsampled to `--record-step-ms`, because CoreNEURON 8.2 does not transfer the
explicit-interval `Vector.record` form used by NetPyNE.

Example matching the 3 s fixed-seed dynamics run except for unsupported LFP:

```bash
uv run -- mpiexec --use-hwthread-cpus --bind-to hwthread \
  --mca btl self,vader,tcp -n 15 python scripts/run_full_quiet_pt5b.py \
  --backend coreneuron --trial 0 --duration-ms 3000 --record-step-ms 0.1 \
  --record-pt5b-all --no-record-lfp --output-dir OUTPUT_DIR
```

Compare saved runs with:

```bash
uv run python scripts/compare_neuron_backends.py REFERENCE_DIR CANDIDATE_DIR
```
