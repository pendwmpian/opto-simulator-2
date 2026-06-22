#!/usr/bin/env python3
"""Zoomed spike raster grouped and colored by model population."""

from __future__ import annotations

import argparse
import json
import math
import pickle
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


LOCAL_POP_ORDER = [
    "IT2",
    "SOM2",
    "PV2",
    "IT4",
    "IT5A",
    "SOM5A",
    "PV5A",
    "IT5B",
    "PT5B",
    "SOM5B",
    "PV5B",
    "IT6",
    "CT6",
    "SOM6",
    "PV6",
]
LONG_POP_ORDER = ["TPO", "TVL", "S1", "S2", "cM1", "M2", "OC"]

POP_COLORS = {
    "IT2": "#f28e8c",
    "IT4": "#e15759",
    "IT5A": "#b03a2e",
    "IT5B": "#7b241c",
    "PT5B": "#1f4aff",
    "IT6": "#641e16",
    "CT6": "#2ca25f",
    "SOM2": "#f8c471",
    "SOM5A": "#f5b041",
    "SOM5B": "#f39c12",
    "SOM6": "#d68910",
    "PV2": "#f7dc6f",
    "PV5A": "#f4d03f",
    "PV5B": "#f1c40f",
    "PV6": "#d4ac0d",
    "TPO": "#7f7f7f",
    "TVL": "#9467bd",
    "S1": "#8c564b",
    "S2": "#e377c2",
    "cM1": "#17becf",
    "M2": "#bcbd22",
    "OC": "#aec7e8",
}


@dataclass(frozen=True)
class PopRange:
    pop: str
    gid_start: int
    gid_stop: int
    active_n: int
    compact_start: int
    compact_stop: int
    group: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path, help="Run directory containing spikes.npz")
    parser.add_argument("--output", type=Path, required=True, help="Output PNG path")
    parser.add_argument("--start-ms", type=float, default=1900.0)
    parser.add_argument("--stop-ms", type=float, default=2050.0)
    parser.add_argument("--bin-ms", type=float, default=2.0)
    parser.add_argument("--baseline-start-ms", type=float, default=1500.0)
    parser.add_argument("--baseline-stop-ms", type=float, default=1850.0)
    parser.add_argument(
        "--density-file",
        type=Path,
        default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/cellDensity.pkl"),
    )
    return parser.parse_args()


def build_pop_ranges(run_dir: Path, density_file: Path) -> list[PopRange]:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    scale = float(manifest["scale"])
    size_x = 300.0
    size_y = 1350.0
    size_z = 300.0
    volume = size_y / 1e3 * size_x / 1e3 / 2.0 * size_z / 1e3 / 2.0 * math.pi
    layer = {
        "2": [0.1, 0.29],
        "24": [0.1, 0.37],
        "4": [0.29, 0.37],
        "5A": [0.37, 0.47],
        "5B": [0.47, 0.8],
        "6": [0.8, 1.0],
    }

    with density_file.open("rb") as f:
        density = pickle.load(f)["density"]

    pop_density = [
        ("IT2", density[("M1", "E")][0], layer["2"], "local E"),
        ("SOM2", density[("M1", "SOM")][5], layer["24"], "local SOM"),
        ("PV2", density[("M1", "PV")][5], layer["24"], "local PV"),
        ("IT4", density[("M1", "E")][1], layer["4"], "local E"),
        ("IT5A", density[("M1", "E")][2], layer["5A"], "local E"),
        ("SOM5A", density[("M1", "SOM")][2], layer["5A"], "local SOM"),
        ("PV5A", density[("M1", "PV")][2], layer["5A"], "local PV"),
        ("IT5B", 0.5 * density[("M1", "E")][3], layer["5B"], "local E"),
        ("PT5B", 0.5 * density[("M1", "E")][3], layer["5B"], "local E"),
        ("SOM5B", density[("M1", "SOM")][3], layer["5B"], "local SOM"),
        ("PV5B", density[("M1", "PV")][3], layer["5B"], "local PV"),
        ("IT6", 0.5 * density[("M1", "E")][4], layer["6"], "local E"),
        ("CT6", 0.5 * density[("M1", "E")][4], layer["6"], "local E"),
        ("SOM6", density[("M1", "SOM")][4], layer["6"], "local SOM"),
        ("PV6", density[("M1", "PV")][4], layer["6"], "local PV"),
    ]

    ranges: list[PopRange] = []
    gid_cursor = 0
    compact_cursor = 0
    for pop, dens, yrange, group in pop_density:
        n = int(float(dens) * volume * (yrange[1] - yrange[0]))
        ranges.append(
            PopRange(pop, gid_cursor, gid_cursor + n - 1, n, compact_cursor, compact_cursor + n - 1, group)
        )
        gid_cursor += n
        compact_cursor += n

    long_declared_n = int(manifest.get("numCellsLong", 1000))
    long_active_n = int(scale * long_declared_n)
    for pop in LONG_POP_ORDER:
        ranges.append(
            PopRange(
                pop,
                gid_cursor,
                gid_cursor + long_active_n - 1,
                long_active_n,
                compact_cursor,
                compact_cursor + long_active_n - 1,
                "long input",
            )
        )
        gid_cursor += long_declared_n
        compact_cursor += long_active_n

    return ranges


def map_spikes_to_pops(spkid: np.ndarray, ranges: list[PopRange]) -> tuple[np.ndarray, np.ndarray]:
    pop_idx = np.full(spkid.shape, -1, dtype=int)
    compact_y = np.full(spkid.shape, np.nan, dtype=float)
    for i, pr in enumerate(ranges):
        mask = (spkid >= pr.gid_start) & (spkid <= pr.gid_stop)
        pop_idx[mask] = i
        compact_y[mask] = pr.compact_start + (spkid[mask] - pr.gid_start)
    return pop_idx, compact_y


def first_event_times(
    spkt: np.ndarray,
    pop_idx: np.ndarray,
    ranges: list[PopRange],
    start_ms: float,
    stop_ms: float,
    baseline_start_ms: float,
    baseline_stop_ms: float,
    bin_ms: float,
) -> list[dict[str, float | str | int]]:
    bins = np.arange(start_ms, stop_ms + bin_ms, bin_ms)
    base_duration_s = max((baseline_stop_ms - baseline_start_ms) / 1000.0, 1e-9)
    rows = []
    for i, pr in enumerate(ranges):
        win_times = spkt[(pop_idx == i) & (spkt >= start_ms) & (spkt < stop_ms)]
        base_count = np.count_nonzero((pop_idx == i) & (spkt >= baseline_start_ms) & (spkt < baseline_stop_ms))
        base_rate = base_count / base_duration_s / pr.active_n
        counts, _ = np.histogram(win_times, bins=bins)
        rates = counts / (bin_ms / 1000.0) / pr.active_n
        threshold = max(base_rate + 3.0, base_rate * 3.0, 1.0 / (bin_ms / 1000.0) / pr.active_n)
        above = np.where(rates >= threshold)[0]
        onset = float(bins[above[0]]) if above.size else np.nan
        peak_i = int(np.argmax(rates)) if rates.size else -1
        peak_time = float(bins[peak_i] + bin_ms / 2.0) if peak_i >= 0 else np.nan
        rows.append(
            {
                "pop": pr.pop,
                "group": pr.group,
                "n": pr.active_n,
                "spikes": int(win_times.size),
                "baseline_hz_cell": float(base_rate),
                "onset_ms": onset,
                "peak_ms": peak_time,
                "peak_hz_cell": float(rates[peak_i]) if peak_i >= 0 else np.nan,
            }
        )
    rows.sort(key=lambda r: (np.inf if np.isnan(float(r["onset_ms"])) else float(r["onset_ms"]), -int(r["spikes"])))
    return rows


def main() -> None:
    args = parse_args()
    manifest = json.loads((args.run_dir / "manifest.json").read_text())
    ranges = build_pop_ranges(args.run_dir, args.density_file)

    with np.load(args.run_dir / "spikes.npz", allow_pickle=False) as spikes:
        spkt = np.asarray(spikes["spkt"], dtype=float)
        spkid = np.asarray(spikes["spkid"], dtype=int)

    pop_idx, compact_y = map_spikes_to_pops(spkid, ranges)
    valid = pop_idx >= 0
    spkt = spkt[valid]
    spkid = spkid[valid]
    pop_idx = pop_idx[valid]
    compact_y = compact_y[valid]

    win = (spkt >= args.start_ms) & (spkt <= args.stop_ms)
    rows = first_event_times(
        spkt,
        pop_idx,
        ranges,
        args.start_ms,
        args.stop_ms,
        args.baseline_start_ms,
        args.baseline_stop_ms,
        args.bin_ms,
    )

    bins = np.arange(args.start_ms, args.stop_ms + args.bin_ms, args.bin_ms)
    centers = bins[:-1] + args.bin_ms / 2.0

    fig = plt.figure(figsize=(14.5, 10.2), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[4.2, 1.8, 1.8], width_ratios=[3.2, 1.25])
    ax_raster = fig.add_subplot(gs[0, :])
    ax_rate = fig.add_subplot(gs[1, 0], sharex=ax_raster)
    ax_heat = fig.add_subplot(gs[2, 0], sharex=ax_raster)
    ax_table = fig.add_subplot(gs[1:, 1])

    for i, pr in enumerate(ranges):
        mask = win & (pop_idx == i)
        if not np.any(mask):
            continue
        ax_raster.scatter(
            spkt[mask],
            compact_y[mask],
            s=4.5 if pr.group == "long input" else 2.2,
            color=POP_COLORS.get(pr.pop, "#444444"),
            alpha=0.75 if pr.group == "long input" else 0.55,
            linewidths=0,
            label=pr.pop,
        )

    for pr in ranges:
        ax_raster.axhline(pr.compact_start - 0.5, color="#e0e0e0", lw=0.45)
        ax_raster.text(
            args.start_ms - (args.stop_ms - args.start_ms) * 0.01,
            (pr.compact_start + pr.compact_stop) / 2.0,
            pr.pop,
            ha="right",
            va="center",
            fontsize=8,
            color=POP_COLORS.get(pr.pop, "#444444"),
        )
    ax_raster.set_xlim(args.start_ms, args.stop_ms)
    ax_raster.set_ylim(-50, ranges[-1].compact_stop + 50)
    ax_raster.set_yticks([])
    ax_raster.set_ylabel("Population-grouped cells")
    ax_raster.set_title(
        f"Zoomed spike raster by population, {args.start_ms:g}-{args.stop_ms:g} ms "
        f"(trial {manifest.get('trial')}, ihGbar={manifest.get('ihGbar')})"
    )
    ax_raster.grid(axis="x", color="#d0d0d0", lw=0.5, alpha=0.7)

    selected = [
        "TVL",
        "TPO",
        "cM1",
        "M2",
        "IT5B",
        "PT5B",
        "PV5B",
        "SOM5B",
        "IT6",
        "CT6",
    ]
    for pop in selected:
        i = next((j for j, pr in enumerate(ranges) if pr.pop == pop), None)
        if i is None:
            continue
        pr = ranges[i]
        counts, _ = np.histogram(spkt[pop_idx == i], bins=bins)
        rate = counts / (args.bin_ms / 1000.0) / pr.active_n
        ax_rate.plot(centers, rate, lw=1.15, color=POP_COLORS.get(pop, "#444444"), label=pop)
    ax_rate.set_ylabel("Hz/cell")
    ax_rate.set_title(f"Per-population event rate ({args.bin_ms:g} ms bins)")
    ax_rate.grid(color="#d0d0d0", lw=0.5, alpha=0.7)
    ax_rate.legend(ncol=5, fontsize=8, frameon=False, loc="upper left")

    heat_pops = [r["pop"] for r in rows if int(r["spikes"]) > 0][:14]
    heat = []
    heat_labels = []
    for pop in heat_pops:
        i = next(j for j, pr in enumerate(ranges) if pr.pop == pop)
        pr = ranges[i]
        counts, _ = np.histogram(spkt[pop_idx == i], bins=bins)
        heat.append(counts / (args.bin_ms / 1000.0) / pr.active_n)
        heat_labels.append(pop)
    if heat:
        heat_arr = np.asarray(heat)
        image = ax_heat.imshow(
            heat_arr,
            aspect="auto",
            interpolation="nearest",
            extent=[args.start_ms, args.stop_ms, len(heat_labels), 0],
            cmap="magma",
        )
        fig.colorbar(image, ax=ax_heat, label="Hz/cell")
        ax_heat.set_yticks(np.arange(len(heat_labels)) + 0.5)
        ax_heat.set_yticklabels(heat_labels)
    ax_heat.set_xlabel("Time (ms)")
    ax_heat.set_title("Populations ordered by first elevated-rate bin")

    ax_table.axis("off")
    table_rows = [r for r in rows if int(r["spikes"]) > 0][:16]
    col_labels = ["pop", "onset", "peak", "spks", "peak Hz/cell"]
    cell_text = [
        [
            str(r["pop"]),
            "n/a" if np.isnan(float(r["onset_ms"])) else f"{float(r['onset_ms']):.1f}",
            f"{float(r['peak_ms']):.1f}",
            str(int(r["spikes"])),
            f"{float(r['peak_hz_cell']):.1f}",
        ]
        for r in table_rows
    ]
    table = ax_table.table(cellText=cell_text, colLabels=col_labels, loc="center", cellLoc="right")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.35)
    for row_i, r in enumerate(table_rows, start=1):
        color = POP_COLORS.get(str(r["pop"]), "#444444")
        table[row_i, 0].get_text().set_color(color)
    ax_table.set_title(
        "Temporal order\nrate onset, not causal proof",
        fontsize=11,
        pad=12,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)

    print(f"wrote {args.output}")
    print("population order by elevated-rate onset:")
    for r in table_rows:
        onset = "n/a" if np.isnan(float(r["onset_ms"])) else f"{float(r['onset_ms']):.1f}"
        print(
            f"{r['pop']}: onset={onset} ms, peak={float(r['peak_ms']):.1f} ms, "
            f"spikes={int(r['spikes'])}, peak_hz_cell={float(r['peak_hz_cell']):.2f}"
        )


if __name__ == "__main__":
    main()
