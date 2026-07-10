#!/usr/bin/env python3
"""Compare saved NEURON and CoreNEURON runs and their wall times."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def array_metrics(reference: np.ndarray, candidate: np.ndarray, atol: float) -> dict[str, Any]:
    result: dict[str, Any] = {
        "reference_shape": list(reference.shape),
        "candidate_shape": list(candidate.shape),
        "same_shape": reference.shape == candidate.shape,
        "bitwise_equal": np.array_equal(reference, candidate, equal_nan=True),
    }
    if reference.shape != candidate.shape:
        result.update({"within_tolerance": False, "max_abs_error": None, "rmse": None})
        return result
    difference = np.abs(reference.astype(np.float64) - candidate.astype(np.float64))
    result.update({
        "within_tolerance": bool(np.allclose(reference, candidate, rtol=0.0, atol=atol, equal_nan=True)),
        "max_abs_error": float(np.nanmax(difference)) if difference.size else 0.0,
        "rmse": float(np.sqrt(np.nanmean(difference * difference))) if difference.size else 0.0,
        "atol": atol,
    })
    return result


def integration_seconds(run_dir: Path, manifest: dict[str, Any]) -> float | None:
    if "simulation_seconds" in manifest:
        return float(manifest["simulation_seconds"])
    resource = run_dir / "resource_rank0.jsonl"
    if not resource.exists():
        return None
    phases = {
        row["phase"]: float(row["wall_time_epoch"])
        for line in resource.read_text(encoding="utf-8").splitlines()
        if (row := json.loads(line))
    }
    if "after_setup_recording" in phases and "after_run_sim" in phases:
        return phases["after_run_sim"] - phases["after_setup_recording"]
    return None


def compare(reference_dir: Path, candidate_dir: Path) -> dict[str, Any]:
    reference_manifest = read_json(reference_dir / "manifest.json")
    candidate_manifest = read_json(candidate_dir / "manifest.json")
    parameter_keys = (
        "duration_ms", "scale", "ihGbar", "seeds", "nhost", "record_step_ms",
        "record_pt5b_all", "target_gid", "pt5b_gid_count", "upstream_commit",
    )
    parameters = {
        key: {
            "reference": reference_manifest.get(key),
            "candidate": candidate_manifest.get(key),
            "equal": reference_manifest.get(key) == candidate_manifest.get(key),
        }
        for key in parameter_keys
    }
    arrays: dict[str, Any] = {}
    for filename, keys, tolerance in (
        ("spikes.npz", ("spkt", "spkid"), 0.0),
        ("target_vm.npz", ("t_ms", "v_mV", "target_gid"), 1e-6),
        ("pt5b_vm.npz", ("t_ms", "gids", "v_mV"), 1e-6),
        ("lfp.npz", ("t_ms", "lfp_mV", "electrode_um", "labels"), 1e-9),
    ):
        reference_path = reference_dir / filename
        candidate_path = candidate_dir / filename
        if not reference_path.exists() or not candidate_path.exists():
            arrays[filename] = {
                "available": False,
                "reference_exists": reference_path.exists(),
                "candidate_exists": candidate_path.exists(),
            }
            continue
        with np.load(reference_path) as reference, np.load(candidate_path) as candidate:
            file_result = {key: array_metrics(reference[key], candidate[key], tolerance) for key in keys}
            arrays[filename] = {"available": True, **file_result}

    reference_total = float(reference_manifest["elapsed_seconds"])
    candidate_total = float(candidate_manifest["elapsed_seconds"])
    reference_sim = integration_seconds(reference_dir, reference_manifest)
    candidate_sim = integration_seconds(candidate_dir, candidate_manifest)
    return {
        "reference_dir": str(reference_dir),
        "candidate_dir": str(candidate_dir),
        "parameters": parameters,
        "all_parameters_equal": all(item["equal"] for item in parameters.values()),
        "arrays": arrays,
        "timing": {
            "reference_total_seconds": reference_total,
            "candidate_total_seconds": candidate_total,
            "total_speedup": reference_total / candidate_total,
            "reference_simulation_seconds": reference_sim,
            "candidate_simulation_seconds": candidate_sim,
            "simulation_speedup": (
                reference_sim / candidate_sim if reference_sim is not None and candidate_sim else None
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("candidate_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(args.reference_dir.resolve(), args.candidate_dir.resolve())
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
