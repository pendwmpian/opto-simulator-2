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
from typing import Any

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "external" / "M1_NetPyNE_CellReports_2023"
SIM_DIR = UPSTREAM / "sim"
DEFAULT_SEEDS = {"conn": 4321, "stim": 1234, "loc": 4321}
LFP_LABELS = ["L5A_600um", "upper_L5B_800um", "lower_L5B_1000um"]
CORENEURON_RANGE_GLOBALS = {
    "tadj_Nca",
    "hinf_cancr", "minf_cancr", "s_inf_cancr",
    "hinf_catcb", "minf_catcb",
    "minf_ch_CavL", "mtau_ch_CavL",
    "ninf_ch_KvAngf", "linf_ch_KvAngf", "taul_ch_KvAngf", "taun_ch_KvAngf",
    "oinf_ch_KvCaB", "otau_ch_KvCaB",
    "minf_naz", "hinf_naz", "mtau_naz", "htau_naz", "tadj_naz",
}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(jsonable(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def effective_config(cfg, *, trial, record_pt5b_all, build_only):
    keys = [
        "dt", "hParams", "scale", "sizeX", "sizeY", "sizeZ", "cellmod", "ihModel",
        "ihGbar", "ihGbarBasal", "ihlkc", "ihlkcBasal", "ihlkcBelowSoma", "ihlke",
        "ihSlope", "KgbarFactor", "EEGain", "EIGain", "IEGain", "IIGain", "IPTGain",
        "IFullGain", "IEweights", "IIweights", "weightNorm", "weightNormThreshold",
        "addConn", "addSubConn", "addLongConn", "numCellsLong", "noiseLong", "delayLong",
        "weightLong", "startLong", "ratesLong", "addPulses", "pulse", "pulse2",
        "addIClamp", "addNetStim", "recordStep", "recordLFP", "saveLFPPops",
        "saveLFPCells", "compactConnFormat",
    ]
    result = {key: getattr(cfg, key, None) for key in keys}
    result.update({
        "duration_ms": float(cfg.duration),
        "trial": int(trial),
        "seeds": dict(cfg.seeds),
        "seed_policy": "fixed_explicit; trial does not modify seeds",
        "dynamic_hd_gbar_modification": False,
        "record_pt5b_all": bool(record_pt5b_all),
        "build_only": bool(build_only),
        "backend": "coreneuron" if bool(getattr(cfg, "coreneuron", False)) else "neuron",
    })
    return jsonable(result)


def require_lfp_compatible_ptrvector(h) -> None:
    """Reject NEURON builds whose PtrVector API leaks during NetPyNE LFP runs."""
    if not hasattr(h.PtrVector(1), "ptr_update_callback"):
        raise RuntimeError(
            "LFP recording requires NEURON PtrVector.ptr_update_callback; "
            "use the project-pinned NEURON 8.2.7 environment"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trial", type=int, default=0, help="Metadata only; does not alter seeds")
    parser.add_argument("--duration-ms", type=float, default=3000.0)
    parser.add_argument(
        "--input-duration-ms",
        type=float,
        help="Generate duration-dependent inputs for this many ms, then simulate --duration-ms",
    )
    parser.add_argument("--record-step-ms", type=float, default=0.1)
    parser.add_argument("--conn-seed", type=int, default=DEFAULT_SEEDS["conn"])
    parser.add_argument("--stim-seed", type=int, default=DEFAULT_SEEDS["stim"])
    parser.add_argument("--loc-seed", type=int, default=DEFAULT_SEEDS["loc"])
    parser.add_argument("--target-gid", type=int, help="Defaults to the lowest PT5B GID")
    parser.add_argument(
        "--diagnostic-gids",
        type=int,
        nargs="+",
        default=[],
        help="Record soma and actual spike-source voltage plus cell metadata for these GIDs",
    )
    parser.add_argument(
        "--disable-hh-simple-mechanisms",
        nargs="+",
        default=[],
        help="Diagnostic only: remove these mechanisms from HH_simple cell rules",
    )
    parser.add_argument("--upstream-dir", type=Path, default=UPSTREAM)
    parser.add_argument("--scale", type=float)
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--record-pt5b-all", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--record-lfp", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--single-cell-pops", action="store_true")
    parser.add_argument("--backend", choices=("coreneuron", "neuron"), default="coreneuron")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if any(args.output_dir.iterdir()):
        raise FileExistsError(f"Output directory must be empty: {args.output_dir}")
    input_duration_ms = args.duration_ms if args.input_duration_ms is None else args.input_duration_ms
    if args.duration_ms <= 0 or args.record_step_ms <= 0:
        raise ValueError("duration and record step must be positive")
    if input_duration_ms < args.duration_ms:
        raise ValueError("input duration must be at least the simulation duration")
    upstream_dir = args.upstream_dir.resolve()
    sim_dir = upstream_dir / "sim"
    if not sim_dir.exists():
        raise FileNotFoundError(f"Upstream sim directory not found: {sim_dir}")
    os.chdir(sim_dir)
    sys.path[:0] = [str(sim_dir), str(ROOT)]

    from neuron import h
    h.nrnmpi_init()
    import neuron
    from netpyne import __version__ as netpyne_version
    from netpyne import sim

    cfg = load_module("upstream_cfg", sim_dir / "cfg.py").cfg
    overrides = load_module("quiet_overrides", ROOT / "configs" / "full_quiet_pt5b.py")
    seeds = {"conn": args.conn_seed, "stim": args.stim_seed, "loc": args.loc_seed}
    cfg = overrides.apply_run_only_overrides(
        cfg,
        duration_ms=input_duration_ms,
        trial=args.trial,
        output_dir=args.output_dir,
        seeds=seeds,
        record_step_ms=args.record_step_ms,
        record_lfp=args.record_lfp and not args.build_only,
    )
    record_lfp = bool(cfg.recordLFP)
    if record_lfp:
        require_lfp_compatible_ptrvector(h)
    if args.scale is not None:
        cfg.scale = args.scale
    cfg.coreneuron = args.backend == "coreneuron"
    cfg.gpu = False
    if cfg.coreneuron:
        from neuron import coreneuron

        # NetPyNE 1.1.1 enables CoreNEURON only after finitialize. Enabling it
        # before setupRecording is required for Vector.record transfer.
        coreneuron.enable = True
        coreneuron.gpu = False
        if record_lfp:
            raise RuntimeError(
                "CoreNEURON cannot execute NetPyNE 1.1.1's per-timestep Python callbacks "
                "for LFP recording; rerun with --no-record-lfp"
            )
    if args.single_cell_pops:
        cfg.scale = 1.0
        cfg.singleCellPops = 1
        cfg.numCellsLong = 1

    # The upstream netParams module imports cfg from __main__.
    globals()["cfg"] = cfg
    net_params = load_module("upstream_netparams", sim_dir / "netParams.py").netParams
    # These writable GLOBALs become per-instance RANGE variables in the
    # CoreNEURON-compatible MOD sources used by both backends. Do not let
    # NetPyNE restore the old HOC globals imported with the cell rules;
    # INITIAL/rates sets the per-instance values.
    for cell_rule in net_params.cellParams.values():
        rule_globals = cell_rule.get("globals", {})
        for name in CORENEURON_RANGE_GLOBALS:
            rule_globals.pop(name, None)
        cell_models = cell_rule.get("conds", {}).get("cellModel", [])
        if isinstance(cell_models, str):
            cell_models = [cell_models]
        if "HH_simple" in cell_models:
            for sec in cell_rule.get("secs", {}).values():
                for mechanism in args.disable_hh_simple_mechanisms:
                    sec.get("mechs", {}).pop(mechanism, None)
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
    if rank == 0:
        write_json(args.output_dir / "effective_config.json", effective_config(
            cfg,
            trial=args.trial,
            record_pt5b_all=args.record_pt5b_all,
            build_only=args.build_only,
        ))
        write_json(args.output_dir / "manifest.json", {
            "condition": "control_quiet",
            "backend": args.backend,
            "status": "running",
            "trial": args.trial,
            "duration_ms": args.duration_ms,
            "input_generation_duration_ms": input_duration_ms,
            "ihGbar": cfg.ihGbar,
            "dynamic_hd_gbar_modification": False,
            "seeds": dict(cfg.seeds),
            "seed_policy": "fixed_explicit; trial does not modify seeds",
            "nhost": nhost,
            "record_step_ms": cfg.recordStep,
            "lfp_electrodes_um": cfg.recordLFP,
            "ptr_update_callback_available": hasattr(h.PtrVector(1), "ptr_update_callback"),
            "coreneuron_enabled": bool(cfg.coreneuron),
            "command": [sys.executable, *sys.argv],
            "effective_config_file": "effective_config.json",
        })
    sim.net.createPops()
    mark("after_create_pops")
    sim.net.createCells()
    mark("after_create_cells")
    # VecStim spike trains depend on cfg.duration and are now fully materialized.
    # Restore the requested integration horizon for short prefix comparisons.
    cfg.duration = float(args.duration_ms)

    diagnostic_gids = sorted(set(args.diagnostic_gids))
    local_diagnostic_cells = {}
    for cell in sim.net.cells:
        if cell.gid not in diagnostic_gids or not getattr(cell, "secs", None):
            continue
        source_name = None
        source_sec = None
        source_loc = 0.5
        for sec_name, sec in cell.secs.items():
            if "spikeGenLoc" in sec:
                source_name, source_sec = sec_name, sec
                source_loc = float(sec["spikeGenLoc"])
                break
        if source_sec is None:
            for sec_name, sec in cell.secs.items():
                if len(sec.get("topol", {})) == 0:
                    source_name, source_sec = sec_name, sec
                    break
        if source_sec is None:
            source_name, source_sec = next(iter(cell.secs.items()))
        soma_name = "soma" if "soma" in cell.secs else source_name
        threshold = float(source_sec.get("threshold", net_params.defaultThreshold))
        local_diagnostic_cells[int(cell.gid)] = {
            "gid": int(cell.gid),
            "rank": rank,
            "pop": cell.tags.get("pop"),
            "cell_type": cell.tags.get("cellType"),
            "cell_model": cell.tags.get("cellModel"),
            "soma_sec": soma_name,
            "source_sec": source_name,
            "source_loc": source_loc,
            "threshold_mV": threshold,
            "cell": cell,
            "source_sec_obj": source_sec,
        }

    gathered_diagnostic_meta = sim.pc.py_alltoall(
        [
            {
                gid: {key: value for key, value in meta.items() if key not in ("cell", "source_sec_obj")}
                for gid, meta in local_diagnostic_cells.items()
            }
            if host == 0
            else None
            for host in range(nhost)
        ]
    )
    if rank == 0:
        diagnostic_meta = {
            int(gid): meta
            for group in gathered_diagnostic_meta
            if group
            for gid, meta in group.items()
        }
        missing_diagnostic_gids = sorted(set(diagnostic_gids) - set(diagnostic_meta))
        if missing_diagnostic_gids:
            raise ValueError(f"Diagnostic GIDs not found: {missing_diagnostic_gids}")
    else:
        diagnostic_meta = {}
    diagnostic_meta = sim.pc.py_broadcast(diagnostic_meta, 0)

    local_pt5b = sorted(int(c.gid) for c in sim.net.cells if c.tags.get("pop") == "PT5B")
    gathered_pt5b = sim.pc.py_alltoall([local_pt5b if host == 0 else None for host in range(nhost)])
    if rank == 0:
        all_pt5b_gids = sorted(gid for group in gathered_pt5b if group for gid in group)
    else:
        all_pt5b_gids = []
    all_pt5b_gids = list(sim.pc.py_broadcast(all_pt5b_gids, 0))
    if not all_pt5b_gids:
        raise RuntimeError("No PT5B cells were created")
    target_gid = all_pt5b_gids[0] if args.target_gid is None else int(args.target_gid)
    if target_gid not in all_pt5b_gids:
        raise ValueError(f"Target GID {target_gid} is not in the PT5B population")
    cfg.recordCells = all_pt5b_gids if args.record_pt5b_all else [target_gid]

    sim.net.connectCells()
    mark("after_connect_cells")
    sim.net.addStims()
    mark("after_add_stims")
    diagnostic_vectors = {}
    gathered_diagnostic_traces = []
    if not args.build_only:
        if cfg.coreneuron:
            # CoreNEURON 8.2 transfers Vector.record data at the integration dt,
            # but not NetPyNE's Vector.record(ptr, recordStep) form. Ask NetPyNE
            # to allocate the vectors, then re-register this runner's soma trace
            # without an explicit interval and downsample after transfer.
            cfg.recordStep = cfg.dt
        sim.setupRecording()
        if cfg.coreneuron:
            for cell in sim.net.cells:
                key = f"cell_{cell.gid}"
                if key in sim.simData["V_soma"]:
                    sim.simData["V_soma"][key].record(cell.secs["soma"]["hObj"](0.5)._ref_v)
        for gid, meta in local_diagnostic_cells.items():
            soma_vector = h.Vector()
            source_vector = h.Vector()
            soma_vector.record(meta["cell"].secs[meta["soma_sec"]]["hObj"](0.5)._ref_v)
            source_vector.record(meta["source_sec_obj"]["hObj"](meta["source_loc"])._ref_v)
            diagnostic_vectors[gid] = (soma_vector, source_vector)
        cfg.recordStep = args.record_step_ms
        mark("after_setup_recording")
        sim.runSim()
        mark("after_run_sim")
        local_diagnostic_traces = {
            gid: {
                "soma_v_mV": list(soma_vector),
                "source_v_mV": list(source_vector),
            }
            for gid, (soma_vector, source_vector) in diagnostic_vectors.items()
        }
        gathered_diagnostic_traces = sim.pc.py_alltoall(
            [local_diagnostic_traces if host == 0 else None for host in range(nhost)]
        )
        include_entries = ["spkt", "spkid", "V_soma"]
        if record_lfp:
            include_entries.append("LFP")
        sim.gatherData(
            gatherLFP=record_lfp,
            gatherDipole=False,
            gatherOnlySimData=True,
            includeSimDataEntries=include_entries,
            analyze=False,
        )
        mark("after_gather_data")

    if rank == 0:
        all_data = getattr(sim, "allSimData", {})
        spkt = np.asarray(all_data.get("spkt", []), dtype=np.float64)
        spkid = np.asarray(all_data.get("spkid", []), dtype=np.int64)
        seed_fields = {
            f"{name}_seed": np.asarray([value], dtype=np.int64) for name, value in cfg.seeds.items()
        }
        trace = all_data.get("V_soma", {})
        voltage = np.asarray(trace.get(f"cell_{target_gid}", []), dtype=np.float32)
        if cfg.coreneuron:
            stride = round(args.record_step_ms / cfg.dt)
            if stride < 1 or not np.isclose(stride * cfg.dt, args.record_step_ms):
                raise ValueError("CoreNEURON record step must be an integer multiple of dt")
            expected_samples = round(cfg.duration / args.record_step_ms)
            voltage = voltage[::stride][:expected_samples]
        pt5b_vm_count = 0
        vm_sample_count = int(voltage.size)
        lfp_sample_count = 0
        output_files = []
        if not args.build_only:
            np.savez_compressed(
                args.output_dir / "spikes.npz", spkt=spkt, spkid=spkid, **seed_fields
            )
            np.savez_compressed(
                args.output_dir / "target_vm.npz",
                t_ms=np.arange(voltage.size, dtype=np.float64) * cfg.recordStep,
                v_mV=voltage,
                target_gid=np.asarray([target_gid]),
                **seed_fields,
            )
            output_files = ["spikes.npz", "target_vm.npz"]
        if diagnostic_gids and not args.build_only:
            diagnostic_trace_by_gid = {
                int(gid): trace_data
                for group in gathered_diagnostic_traces
                if group
                for gid, trace_data in group.items()
            }
            stride = round(args.record_step_ms / cfg.dt)
            expected_samples = round(cfg.duration / args.record_step_ms)
            soma_rows = []
            source_rows = []
            for gid in diagnostic_gids:
                trace_data = diagnostic_trace_by_gid[gid]
                soma_rows.append(np.asarray(trace_data["soma_v_mV"], dtype=np.float32)[::stride][:expected_samples])
                source_rows.append(np.asarray(trace_data["source_v_mV"], dtype=np.float32)[::stride][:expected_samples])
            diagnostic_samples = min(
                [row.size for row in soma_rows + source_rows], default=0
            )
            soma_vm = np.stack([row[:diagnostic_samples] for row in soma_rows])
            source_vm = np.stack([row[:diagnostic_samples] for row in source_rows])
            np.savez_compressed(
                args.output_dir / "diagnostic_vm.npz",
                t_ms=np.arange(diagnostic_samples, dtype=np.float64) * cfg.recordStep,
                gids=np.asarray(diagnostic_gids, dtype=np.int64),
                soma_v_mV=soma_vm,
                source_v_mV=source_vm,
                ranks=np.asarray([diagnostic_meta[gid]["rank"] for gid in diagnostic_gids], dtype=np.int64),
                pops=np.asarray([diagnostic_meta[gid]["pop"] for gid in diagnostic_gids]),
                cell_types=np.asarray([diagnostic_meta[gid]["cell_type"] for gid in diagnostic_gids]),
                cell_models=np.asarray([diagnostic_meta[gid]["cell_model"] for gid in diagnostic_gids]),
                soma_secs=np.asarray([diagnostic_meta[gid]["soma_sec"] for gid in diagnostic_gids]),
                source_secs=np.asarray([diagnostic_meta[gid]["source_sec"] for gid in diagnostic_gids]),
                source_locs=np.asarray([diagnostic_meta[gid]["source_loc"] for gid in diagnostic_gids]),
                thresholds_mV=np.asarray([diagnostic_meta[gid]["threshold_mV"] for gid in diagnostic_gids]),
                **seed_fields,
            )
            output_files.append("diagnostic_vm.npz")
        if args.record_pt5b_all and not args.build_only:
            pt5b_rows = [np.asarray(trace.get(f"cell_{gid}", []), dtype=np.float32) for gid in all_pt5b_gids]
            if cfg.coreneuron:
                pt5b_rows = [row[::stride][:expected_samples] for row in pt5b_rows]
            max_len = max((row.size for row in pt5b_rows), default=0)
            pt5b_vm = np.full((len(pt5b_rows), max_len), np.nan, dtype=np.float32)
            for row_index, row in enumerate(pt5b_rows):
                pt5b_vm[row_index, :row.size] = row
            np.savez_compressed(
                args.output_dir / "pt5b_vm.npz",
                t_ms=np.arange(max_len, dtype=np.float64) * cfg.recordStep,
                gids=np.asarray(all_pt5b_gids, dtype=np.int64),
                v_mV=pt5b_vm,
                **seed_fields,
            )
            pt5b_vm_count = len(pt5b_rows)
            vm_sample_count = max_len
            output_files.append("pt5b_vm.npz")
        if record_lfp and not args.build_only:
            lfp = np.asarray(all_data.get("LFP", []), dtype=np.float64)
            if lfp.ndim != 2 or lfp.shape[1] != len(cfg.recordLFP):
                raise RuntimeError(f"Unexpected LFP shape: {lfp.shape}")
            np.savez_compressed(
                args.output_dir / "lfp.npz",
                t_ms=np.arange(lfp.shape[0], dtype=np.float64) * cfg.recordStep,
                lfp_mV=lfp,
                electrode_um=np.asarray(cfg.recordLFP, dtype=np.float64),
                labels=np.asarray(LFP_LABELS),
                **seed_fields,
            )
            lfp_sample_count = int(lfp.shape[0])
            output_files.append("lfp.npz")
        commit = subprocess.check_output(
            ["git", "-C", str(upstream_dir), "rev-parse", "HEAD"], text=True
        ).strip()
        manifest = {
            "condition": "control_quiet", "status": "success", "trial": args.trial,
            "backend": args.backend,
            "duration_ms": cfg.duration, "scale": cfg.scale,
            "input_generation_duration_ms": input_duration_ms,
            "ihGbar": getattr(cfg, "ihGbar", None),
            "dynamic_hd_gbar_modification": False,
            "build_only": args.build_only, "seeds": dict(cfg.seeds),
            "seed_policy": "fixed_explicit; trial does not modify seeds",
            "target_gid": target_gid, "pt5b_gid_count": len(all_pt5b_gids),
            "record_pt5b_all": args.record_pt5b_all, "nhost": nhost,
            "diagnostic_gids": diagnostic_gids,
            "disabled_hh_simple_mechanisms": sorted(set(args.disable_hh_simple_mechanisms)),
            "all_spike_count": int(spkt.size),
            "recorded_pt5b_vm_count": pt5b_vm_count,
            "vm_sample_count": vm_sample_count,
            "lfp_sample_count": lfp_sample_count,
            "elapsed_seconds": time.time() - started,
            "simulation_seconds": float(sim.timingData.get("runTime", 0.0)),
            "upstream_commit": commit, "python": sys.version,
            "platform": platform.platform(), "neuron_version": neuron.__version__,
            "netpyne_version": netpyne_version, "record_step_ms": cfg.recordStep,
            "command": [sys.executable, *sys.argv],
            "lfp_electrodes_um": cfg.recordLFP,
            "ptr_update_callback_available": hasattr(h.PtrVector(1), "ptr_update_callback"),
            "coreneuron_enabled": bool(cfg.coreneuron),
            "effective_config_file": "effective_config.json",
            "output": output_files,
        }
        write_json(args.output_dir / "effective_config.json", effective_config(
            cfg,
            trial=args.trial,
            record_pt5b_all=args.record_pt5b_all,
            build_only=args.build_only,
        ))
        write_json(args.output_dir / "manifest.json", manifest)
        print(json.dumps(manifest, indent=2), flush=True)
    sim.pc.barrier()
    sim.pc.done()
    from mpi4py import MPI
    if not MPI.Is_finalized():
        MPI.Finalize()


if __name__ == "__main__":
    main()
