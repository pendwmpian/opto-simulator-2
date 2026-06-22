#!/usr/bin/env python3
"""Plot all recorded spikes as a raster, optionally highlighting PT5B cells."""

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
    parser.add_argument("run_dir", type=Path, help="Run directory containing spikes.npz")
    parser.add_argument("--output", type=Path, required=True, help="Output PNG path")
    parser.add_argument("--bin-ms", type=float, default=25.0, help="Population-rate bin width")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spikes_path = args.run_dir / "spikes.npz"
    manifest_path = args.run_dir / "manifest.json"
    pt5b_path = args.run_dir / "pt5b_vm.npz"

    with np.load(spikes_path, allow_pickle=False) as spikes:
        spkt = np.asarray(spikes["spkt"], dtype=float)
        spkid = np.asarray(spikes["spkid"], dtype=int)

    if spkt.size == 0:
        raise ValueError(f"No spikes found in {spikes_path}")

    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    pt5b_gids = np.array([], dtype=int)
    if pt5b_path.exists():
        with np.load(pt5b_path, allow_pickle=False) as pt5b:
            pt5b_gids = np.asarray(pt5b["gids"], dtype=int)
    pt5b_mask = np.isin(spkid, pt5b_gids) if pt5b_gids.size else np.zeros(spkid.shape, dtype=bool)

    duration_ms = float(manifest.get("duration_ms", np.nanmax(spkt)))
    bins = np.arange(0.0, duration_ms + args.bin_ms, args.bin_ms)
    all_counts, edges = np.histogram(spkt, bins=bins)
    pt5b_counts, _ = np.histogram(spkt[pt5b_mask], bins=bins)
    centers = edges[:-1] + np.diff(edges) / 2.0
    all_rate_hz = all_counts / (args.bin_ms / 1000.0)
    pt5b_rate_hz = pt5b_counts / (args.bin_ms / 1000.0)

    fig = plt.figure(figsize=(12.5, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[4.2, 1.2])
    ax_raster = fig.add_subplot(gs[0])
    ax_rate = fig.add_subplot(gs[1], sharex=ax_raster)

    ax_raster.scatter(spkt[~pt5b_mask], spkid[~pt5b_mask], s=0.8, c="#303030", alpha=0.18, linewidths=0)
    if np.any(pt5b_mask):
        ax_raster.scatter(spkt[pt5b_mask], spkid[pt5b_mask], s=2.5, c="#c62828", alpha=0.75, linewidths=0)
        ax_raster.axhspan(int(pt5b_gids.min()), int(pt5b_gids.max()), color="#c62828", alpha=0.05, linewidth=0)

    title = (
        f"All recorded spikes, trial {manifest.get('trial', 'unknown')}, "
        f"ihGbar={manifest.get('ihGbar', 'unknown')}"
    )
    ax_raster.set_title(title)
    ax_raster.set_ylabel("GID")
    ax_raster.set_xlim(0.0, duration_ms)
    ax_raster.grid(axis="x", color="#d0d0d0", linewidth=0.5, alpha=0.7)

    ax_rate.plot(centers, all_rate_hz, color="#303030", lw=1.0, label="all spike events / s")
    if np.any(pt5b_mask):
        ax_rate.plot(centers, pt5b_rate_hz, color="#c62828", lw=1.0, label="PT5B spike events / s")
    ax_rate.set_xlabel("Time (ms)")
    ax_rate.set_ylabel("Events/s")
    ax_rate.grid(color="#d0d0d0", linewidth=0.5, alpha=0.7)
    ax_rate.legend(loc="upper right", frameon=False)

    unique_spiking_gids = np.unique(spkid).size
    summary = [
        f"spikes={spkt.size:,}",
        f"spiking gids={unique_spiking_gids:,}",
        f"gid range={int(spkid.min())}-{int(spkid.max())}",
    ]
    if pt5b_gids.size:
        summary.append(f"PT5B spikes={int(pt5b_mask.sum()):,} / cells={pt5b_gids.size:,}")
    ax_raster.text(
        0.01,
        0.98,
        "\n".join(summary),
        transform=ax_raster.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.86, "edgecolor": "#bbbbbb"},
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)
    print(f"wrote {args.output}")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
