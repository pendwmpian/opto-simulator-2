#!/usr/bin/env python3
"""Create a compact overview plot for all recorded PT5B soma Vm traces."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path, help="Run directory containing pt5b_vm.npz")
    parser.add_argument("--output", type=Path, required=True, help="Output PNG path")
    parser.add_argument("--analysis-start-ms", type=float, default=1000.0)
    parser.add_argument("--max-heatmap-time-points", type=int, default=1200)
    return parser.parse_args()


def spike_counts_for_gids(run_dir: Path, gids: np.ndarray) -> np.ndarray:
    spikes_path = run_dir / "spikes.npz"
    if not spikes_path.exists():
        return np.zeros(gids.shape, dtype=int)
    with np.load(spikes_path, allow_pickle=False) as spikes:
        counts = Counter(spikes["spkid"].astype(int))
    return np.asarray([counts.get(int(gid), 0) for gid in gids], dtype=int)


def decimate_time(t_ms: np.ndarray, v_mV: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    if t_ms.size <= max_points:
        return t_ms, v_mV
    step = int(np.ceil(t_ms.size / max_points))
    return t_ms[::step], v_mV[:, ::step]


def main() -> None:
    args = parse_args()
    vm_path = args.run_dir / "pt5b_vm.npz"
    manifest_path = args.run_dir / "manifest.json"

    with np.load(vm_path, allow_pickle=False) as data:
        t_ms = np.asarray(data["t_ms"], dtype=float)
        gids = np.asarray(data["gids"], dtype=int)
        v_mV = np.asarray(data["v_mV"], dtype=float)

    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    analysis_mask = t_ms >= args.analysis_start_ms
    if not np.any(analysis_mask):
        raise ValueError("analysis start is later than the recorded trace")

    post_v = v_mV[:, analysis_mask]
    mean_v = post_v.mean(axis=1)
    sd_v = post_v.std(axis=1)
    max_v = v_mV.max(axis=1)
    min_v = v_mV.min(axis=1)
    spike_counts = spike_counts_for_gids(args.run_dir, gids)

    sort_idx = np.lexsort((-spike_counts, mean_v))
    t_heat, v_heat = decimate_time(t_ms, v_mV[sort_idx], args.max_heatmap_time_points)

    representative_indices = []
    target_gid = manifest.get("target_gid")
    if target_gid in set(gids.tolist()):
        representative_indices.append(int(np.where(gids == int(target_gid))[0][0]))
    for idx in [int(np.argmax(spike_counts)), int(np.argmax(mean_v)), int(np.argmin(mean_v))]:
        if idx not in representative_indices:
            representative_indices.append(idx)

    fig = plt.figure(figsize=(12.5, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[2.5, 1.0, 1.0])
    ax_heat = fig.add_subplot(gs[0, :])
    ax_trace = fig.add_subplot(gs[1, :])
    ax_hist = fig.add_subplot(gs[2, 0])
    ax_scatter = fig.add_subplot(gs[2, 1])

    image = ax_heat.imshow(
        v_heat,
        aspect="auto",
        interpolation="nearest",
        extent=[t_heat[0], t_heat[-1], 0, len(gids)],
        origin="lower",
        cmap="viridis",
        vmin=-80,
        vmax=-45,
    )
    fig.colorbar(image, ax=ax_heat, label="Vm (mV)")
    ax_heat.set_title(
        f"All recorded PT5B soma Vm, trial {manifest.get('trial', 'unknown')}, "
        f"ihGbar={manifest.get('ihGbar', 'unknown')}"
    )
    ax_heat.set_xlabel("Time (ms)")
    ax_heat.set_ylabel("PT5B cells sorted by post-1s mean Vm")

    colors = ["#1769aa", "#c62828", "#2e7d32", "#6a1b9a"]
    for color, idx in zip(colors, representative_indices):
        label = (
            f"GID {int(gids[idx])}: spikes={int(spike_counts[idx])}, "
            f"mean={mean_v[idx]:.1f} mV"
        )
        ax_trace.plot(t_ms, v_mV[idx], lw=0.8, color=color, label=label)
    ax_trace.axvspan(t_ms[0], min(args.analysis_start_ms, t_ms[-1]), color="#9e9e9e", alpha=0.14, lw=0)
    ax_trace.set_xlabel("Time (ms)")
    ax_trace.set_ylabel("Vm (mV)")
    ax_trace.legend(loc="upper right", frameon=False, fontsize=8)
    ax_trace.grid(axis="y", color="#d0d0d0", lw=0.5, alpha=0.7)

    ax_hist.hist(mean_v, bins=40, color="#1769aa", alpha=0.8)
    ax_hist.set_xlabel(f"Mean Vm after {args.analysis_start_ms:g} ms (mV)")
    ax_hist.set_ylabel("Cells")

    ax_scatter.scatter(mean_v, spike_counts, s=10, color="#333333", alpha=0.65)
    ax_scatter.set_xlabel(f"Mean Vm after {args.analysis_start_ms:g} ms (mV)")
    ax_scatter.set_ylabel("Spike count")
    ax_scatter.grid(color="#d0d0d0", lw=0.5, alpha=0.6)

    summary = (
        f"cells={len(gids)}, duration={t_ms[-1] - t_ms[0] + np.median(np.diff(t_ms)):.1f} ms\n"
        f"mean Vm range={mean_v.min():.2f} to {mean_v.max():.2f} mV\n"
        f"spike count range={spike_counts.min()} to {spike_counts.max()}\n"
        f"Vm full range={min_v.min():.2f} to {max_v.max():.2f} mV"
    )
    fig.text(0.01, 0.01, summary, fontsize=9, va="bottom", ha="left")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
