#!/usr/bin/env python3
"""Run the full Dura-Bernal M1 network with minimal retained output."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "external" / "M1_NetPyNE_CellReports_2023"
SIM_DIR = UPSTREAM / "sim"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trial", type=int, default=0, choices=range(3))
    parser.add_argument("--duration-ms", type=float, default=1000.0)
    parser.add_argument("--scale", type=float)
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--record-pt5b-all", action="store_true")
    parser.add_argument("--single-cell-pops", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(SIM_DIR)
    sys.path[:0] = [str(SIM_DIR), str(ROOT)]

    from neuron import h
    h.nrnmpi_init()
    import neuron
    from netpyne import __version__ as netpyne_version
    from netpyne import sim

    cfg = load_module("upstream_cfg", SIM_DIR / "cfg.py").cfg
    overrides = load_module("quiet_overrides", ROOT / "configs" / "full_quiet_pt5b.py")
    cfg = overrides.apply_run_only_overrides(
        cfg, duration_ms=args.duration_ms, trial=args.trial, output_dir=args.output_dir
    )
    if args.scale is not None:
        cfg.scale = args.scale
    if args.single_cell_pops:
        cfg.scale = 1.0
        cfg.singleCellPops = 1
        cfg.numCellsLong = 1

    # The upstream netParams module imports cfg from __main__.
    globals()["cfg"] = cfg
    net_params = load_module("upstream_netparams", SIM_DIR / "netParams.py").netParams
    rank = 0
    nhost = 1
    resource_path = args.output_dir / f"resource_rank{rank}.jsonl"

    def mark(phase):
        process = psutil.Process()
        row = {
            "phase": phase,
            "wall_time_epoch": time.time(),
            "rank": rank,
            "rss_bytes": process.memory_info().rss,
            "system_available_bytes": psutil.virtual_memory().available,
            "neuron_t_ms": float(h.t),
        }
        with resource_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(row) + "\n")
        print(
            f"RESOURCE phase={phase} rank={rank} rss_gib={row['rss_bytes']/2**30:.3f} "
            f"avail_gib={row['system_available_bytes']/2**30:.3f}", flush=True
        )
        if hasattr(sim, "pc"):
            sim.pc.barrier()

    started = time.time()
    mark("before_initialize")
    sim.initialize(simConfig=cfg, netParams=net_params)
    rank = int(sim.rank)
    nhost = int(sim.nhosts)
    resource_path = args.output_dir / f"resource_rank{rank}.jsonl"
    sim.pc.timeout(0)
    mark("after_initialize")
    sim.net.createPops()
    mark("after_create_pops")
    sim.net.createCells()
    mark("after_create_cells")

    local_pt5b = sorted(int(c.gid) for c in sim.net.cells if c.tags.get("pop") == "PT5B")
    gathered_pt5b = sim.pc.py_alltoall([local_pt5b] * nhost)
    all_pt5b_gids = sorted(gid for group in gathered_pt5b for gid in group)
    target_gid = all_pt5b_gids[0]
    cfg.recordCells = all_pt5b_gids if args.record_pt5b_all else [target_gid]

    sim.net.connectCells()
    mark("after_connect_cells")
    sim.net.addStims()
    mark("after_add_stims")
    sim.setupRecording()
    mark("after_setup_recording")
    if not args.build_only:
        sim.runSim()
        mark("after_run_sim")
        sim.gatherData(gatherLFP=False, gatherDipole=False)
        mark("after_gather_data")

    if rank == 0:
        all_data = getattr(sim, "allSimData", {})
        spkt = np.asarray(all_data.get("spkt", []), dtype=np.float64)
        spkid = np.asarray(all_data.get("spkid", []), dtype=np.int64)
        np.savez_compressed(args.output_dir / "spikes.npz", spkt=spkt, spkid=spkid)
        trace = all_data.get("V_soma", {})
        voltage = np.asarray(trace.get(f"cell_{target_gid}", []), dtype=np.float32)
        t_ms = np.arange(voltage.size, dtype=np.float64) * cfg.recordStep
        np.savez_compressed(
            args.output_dir / "target_vm.npz",
            t_ms=t_ms,
            v_mV=voltage,
            target_gid=np.asarray([target_gid]),
        )
        if args.record_pt5b_all:
            pt5b_rows = [np.asarray(trace.get(f"cell_{gid}", []), dtype=np.float32) for gid in all_pt5b_gids]
            max_len = max((row.size for row in pt5b_rows), default=0)
            pt5b_vm = np.full((len(pt5b_rows), max_len), np.nan, dtype=np.float32)
            for row_index, row in enumerate(pt5b_rows):
                pt5b_vm[row_index, :row.size] = row
            np.savez_compressed(
                args.output_dir / "pt5b_vm.npz",
                t_ms=np.arange(max_len, dtype=np.float64) * cfg.recordStep,
                gids=np.asarray(all_pt5b_gids, dtype=np.int64),
                v_mV=pt5b_vm,
            )
        commit = subprocess.check_output(
            ["git", "-C", str(UPSTREAM), "rev-parse", "HEAD"], text=True
        ).strip()
        manifest = {
            "condition": "control_quiet", "trial": args.trial,
            "duration_ms": cfg.duration, "scale": cfg.scale,
            "ihGbar": getattr(cfg, "ihGbar", None),
            "build_only": args.build_only, "seeds": dict(cfg.seeds),
            "target_gid": target_gid, "pt5b_gid_count": len(all_pt5b_gids),
            "record_pt5b_all": args.record_pt5b_all, "nhost": nhost,
            "elapsed_seconds": time.time() - started,
            "upstream_commit": commit, "python": sys.version,
            "platform": platform.platform(), "neuron_version": neuron.__version__,
            "netpyne_version": netpyne_version, "record_step_ms": cfg.recordStep,
            "output": ["all_spikes", "target_pt5b_soma_vm"] + (["all_pt5b_soma_vm"] if args.record_pt5b_all else []),
        }
        (args.output_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(manifest, indent=2), flush=True)
    sim.pc.barrier()
    sim.pc.done()
    from mpi4py import MPI
    if not MPI.Is_finalized():
        MPI.Finalize()


if __name__ == "__main__":
    main()
