#!/usr/bin/env python3
"""Plot the recorded target-cell soma membrane potential from a full-model run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path, help="Run directory containing target_vm.npz")
    parser.add_argument("--output", type=Path, required=True, help="Output PNG or PDF path")
    parser.add_argument(
        "--analysis-start-ms",
        type=float,
        default=100.0,
        help="Start of the summary-statistics interval (default: 100 ms)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    trace_path = args.run_dir / "target_vm.npz"
    manifest_path = args.run_dir / "manifest.json"
    spikes_path = args.run_dir / "spikes.npz"

    with np.load(trace_path, allow_pickle=False) as data:
        t_ms = np.asarray(data["t_ms"], dtype=float)
        v_mv = np.asarray(data["v_mV"], dtype=float)
        target_gid = int(np.asarray(data["target_gid"]).reshape(-1)[0])

    if t_ms.size == 0 or t_ms.shape != v_mv.shape:
        raise ValueError(f"Invalid or empty Vm trace in {trace_path}")
    if not np.all(np.isfinite(t_ms)) or not np.all(np.isfinite(v_mv)):
        raise ValueError(f"Non-finite values found in {trace_path}")

    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    target_spikes = None
    if spikes_path.exists():
        with np.load(spikes_path, allow_pickle=False) as spikes:
            target_spikes = int(np.count_nonzero(spikes["spkid"] == target_gid))
    analysis_mask = t_ms >= args.analysis_start_ms
    if not np.any(analysis_mask):
        raise ValueError("analysis start is later than the recorded trace")

    analysis_v = v_mv[analysis_mask]
    mean_mv = float(np.mean(analysis_v))
    sd_mv = float(np.std(analysis_v))
    minimum_mv = float(np.min(v_mv))
    maximum_mv = float(np.max(v_mv))

    fig, ax = plt.subplots(figsize=(10.5, 4.8), constrained_layout=True)
    ax.plot(t_ms, v_mv, color="#1769aa", linewidth=0.9)
    if args.analysis_start_ms > t_ms[0]:
        ax.axvspan(
            t_ms[0],
            min(args.analysis_start_ms, t_ms[-1]),
            color="#9e9e9e",
            alpha=0.18,
            linewidth=0,
            label="Initialization interval",
        )
    ax.hlines(
        mean_mv,
        max(args.analysis_start_ms, t_ms[0]),
        t_ms[-1],
        color="#c62828",
        linestyle="--",
        linewidth=1.2,
        label=f"Mean after {args.analysis_start_ms:g} ms: {mean_mv:.2f} mV",
    )

    condition = manifest.get("condition", "recorded condition")
    trial = manifest.get("trial", "unknown")
    ax.set_title(f"Recorded PT5B soma membrane potential (GID {target_gid})")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Membrane potential (mV)")
    ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
    ax.legend(loc="lower right", frameon=False)

    summary_lines = [
        f"condition: {condition}, trial: {trial}",
        f"post-{args.analysis_start_ms:g} ms SD: {sd_mv:.2f} mV",
        f"full-trace range: {minimum_mv:.2f} to {maximum_mv:.2f} mV",
    ]
    if target_spikes is not None:
        summary_lines.append(f"target spikes: {target_spikes}")
    summary = "\n".join(summary_lines)
    ax.text(
        0.015,
        0.97,
        summary,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "alpha": 0.88, "edgecolor": "#bbbbbb"},
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)

    print(f"wrote {args.output}")
    print(f"target_gid={target_gid}")
    print(f"post_start_mean_mV={mean_mv:.6f}")
    print(f"post_start_sd_mV={sd_mv:.6f}")
    print(f"full_range_mV={minimum_mv:.6f},{maximum_mv:.6f}")


if __name__ == "__main__":
    main()
