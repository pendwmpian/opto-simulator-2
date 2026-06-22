# Full quiet run: fixed seeds and L5 LFP

The full-network quiet runner is `scripts/run_full_quiet_pt5b.py`.

## Fixed scientific condition

- Duration: 3000 ms by default.
- PT `ihGbar`: 0.75 for the entire run.
- Dynamic `hd.gbar` modification: disabled.
- Explicit pulse, IClamp, and test NetStim: disabled.
- Default seeds are fixed and are not derived from the trial number:
  - connection: 4321
  - stimulation/background input: 1234
  - location: 4321

To compare conditions, reuse the same three seed arguments. A different
`--trial` value only changes the metadata label.

## LFP recording

Three electrodes reproduce the L5 locations used by the upstream analysis:

- `[150, 600, 150]` um: L5A
- `[150, 800, 150]` um: upper L5B
- `[150, 1000, 150]` um: lower L5B

Only the summed LFP is saved. Per-cell and per-population LFP contributions are
disabled. MPI rank-local LFP arrays are summed during gather.

## Example command

```bash
uv run -- mpiexec --use-hwthread-cpus --bind-to hwthread --mca btl self,vader,tcp -n 16 \
  python scripts/run_full_quiet_pt5b.py \
  --trial 0 \
  --duration-ms 3000 \
  --conn-seed 4321 \
  --stim-seed 1234 \
  --loc-seed 4321 \
  --output-dir docs-temp/full-model-runs/full-quiet-3s-fixed-seed-lfp-trial0
```

## Output

- `manifest.json`: status, fixed seeds, versions, command, and output summary.
- `effective_config.json`: minimal comparison-relevant effective configuration.
- `spikes.npz`: all gathered spike times and GIDs.
- `target_vm.npz`: target PT5B soma Vm.
- `pt5b_vm.npz`: all PT5B soma Vm unless `--no-record-pt5b-all` is used.
- `lfp.npz`: time, summed LFP, electrode coordinates, and depth labels.
- `resource_rank*.jsonl`: per-rank phase and memory markers.

The three seed values are also embedded in every `.npz` output as
`conn_seed`, `stim_seed`, and `loc_seed`.

The runner does not serialize the full NetParams, cell sections, connections,
or per-cell/per-population LFP data.
