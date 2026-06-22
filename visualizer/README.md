# Full-network result visualizers

These scripts visualize the compact outputs produced by
`scripts/run_full_quiet_pt5b.py`.

## Interactive PT5B Vm viewer

```bash
uv run python visualizer/serve_pt5b_vm_viewer.py RUN_DIR
```

## Static plots

```bash
uv run python visualizer/plot_recorded_target_vm.py RUN_DIR --output target_vm.png
uv run python visualizer/plot_pt5b_vm_population.py RUN_DIR --output pt5b_population.png
uv run python visualizer/plot_spike_raster.py RUN_DIR --output spike_raster.png
uv run python visualizer/plot_spike_raster_by_population_zoom.py RUN_DIR \
  --start-ms 1800 --stop-ms 2200 --output spike_raster_zoom.png
```

The population-grouped raster uses the upstream model population metadata and
reports temporal rate onset as descriptive ordering, not proof of causality.
