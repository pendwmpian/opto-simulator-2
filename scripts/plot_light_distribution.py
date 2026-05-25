from __future__ import annotations

import argparse
import csv
import math
import os
import pickle
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def load_pickle(path: Path):
    with path.open("rb") as f:
        try:
            return pickle.load(f)
        except UnicodeDecodeError:
            f.seek(0)
            return pickle.load(f, encoding="latin1")


def section_kind(name: str) -> str:
    if name.startswith("soma"):
        return "soma"
    if name.startswith("axon"):
        return "axon"
    if name.startswith("apic") or name.startswith("Adend"):
        return "apical_dendrite"
    if name.startswith("dend") or name.startswith("Bdend"):
        return "basal_dendrite"
    return "other"


def collect_segments(rule):
    rows = []
    surface_y = max(
        float(point[1])
        for sec in rule["secs"].values()
        for point in sec.get("geom", {}).get("pt3d", [])
    )
    for sec_name, sec in rule["secs"].items():
        kind = section_kind(sec_name)
        pts = sec.get("geom", {}).get("pt3d", [])
        if len(pts) < 2:
            continue
        for p0, p1 in zip(pts[:-1], pts[1:]):
            x0, y0, z0, d0 = [float(v) for v in p0]
            x1, y1, z1, d1 = [float(v) for v in p1]
            length_um = math.dist((x0, y0, z0), (x1, y1, z1))
            if length_um == 0:
                continue
            diam_um = max(0.01, 0.5 * (d0 + d1))
            area_um2 = math.pi * diam_um * length_um
            mx = 0.5 * (x0 + x1)
            my = 0.5 * (y0 + y1)
            mz = 0.5 * (z0 + z1)
            rows.append(
                {
                    "sec": sec_name,
                    "kind": kind,
                    "p0": (x0, y0, z0),
                    "p1": (x1, y1, z1),
                    "mid": (mx, my, mz),
                    "depth_um": max(0.0, surface_y - my),
                    "area_um2": area_um2,
                }
            )
    return rows


def add_light(rows, surface_irradiance_mw_mm2: float, mu_eff_mm_inv: float):
    for row in rows:
        depth_mm = row["depth_um"] / 1000.0
        transmission = math.exp(-mu_eff_mm_inv * depth_mm)
        local_irradiance = surface_irradiance_mw_mm2 * transmission
        area_mm2 = row["area_um2"] * 1e-6
        row["transmission"] = transmission
        row["irradiance_mw_mm2"] = local_irradiance
        row["intercepted_power_mw"] = local_irradiance * area_mm2


def summarize(rows):
    groups: dict[str, dict[str, float]] = {}
    for row in rows:
        if row["kind"] == "axon":
            continue
        g = groups.setdefault(
            row["kind"],
            {
                "area_um2": 0.0,
                "intercepted_power_mw": 0.0,
                "weighted_irradiance_sum": 0.0,
                "min_irradiance_mw_mm2": float("inf"),
                "max_irradiance_mw_mm2": 0.0,
            },
        )
        area = row["area_um2"]
        irr = row["irradiance_mw_mm2"]
        g["area_um2"] += area
        g["intercepted_power_mw"] += row["intercepted_power_mw"]
        g["weighted_irradiance_sum"] += irr * area
        g["min_irradiance_mw_mm2"] = min(g["min_irradiance_mw_mm2"], irr)
        g["max_irradiance_mw_mm2"] = max(g["max_irradiance_mw_mm2"], irr)
    total_power = sum(g["intercepted_power_mw"] for g in groups.values())
    summary = []
    for kind, g in sorted(groups.items()):
        area = g["area_um2"]
        power = g["intercepted_power_mw"]
        summary.append(
            {
                "part": kind,
                "area_um2": area,
                "area_weighted_irradiance_mw_mm2": g["weighted_irradiance_sum"] / area,
                "min_irradiance_mw_mm2": g["min_irradiance_mw_mm2"],
                "max_irradiance_mw_mm2": g["max_irradiance_mw_mm2"],
                "intercepted_power_mw": power,
                "fraction_of_non_axon_power": power / total_power if total_power else 0.0,
            }
        )
    return summary


def set_equal_3d(ax, rows):
    pts = [p for row in rows for p in (row["p0"], row["p1"])]
    xs, ys, zs = zip(*pts)
    cx, cy, cz = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2
    radius = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)) / 2
    ax.set_xlim(cx - radius, cx + radius)
    ax.set_ylim(cy - radius, cy + radius)
    ax.set_zlim(cz - radius, cz + radius)


def plot(rows, out: Path, title: str):
    non_axon = [r for r in rows if r["kind"] != "axon"]
    seg3d = [[r["p0"], r["p1"]] for r in rows]
    vals = np.array([r["irradiance_mw_mm2"] for r in rows])
    norm = matplotlib.colors.LogNorm(vmin=max(vals[vals > 0].min(), 1e-4), vmax=vals.max())
    cmap = "viridis"

    fig = plt.figure(figsize=(14, 6), constrained_layout=True)
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig.add_subplot(1, 2, 2)

    lc3d = Line3DCollection(seg3d, cmap=cmap, norm=norm, linewidths=1.2)
    lc3d.set_array(vals)
    ax3d.add_collection3d(lc3d)
    set_equal_3d(ax3d, rows)
    ax3d.set_xlabel("x (um)")
    ax3d.set_ylabel("y (um; + toward pia locally)")
    ax3d.set_zlabel("z (um)")
    ax3d.view_init(elev=14, azim=-68)
    ax3d.set_title("3D morphology colored by local irradiance")

    seg2d = [[(r["p0"][0], r["p0"][1]), (r["p1"][0], r["p1"][1])] for r in rows]
    lc2d = LineCollection(seg2d, cmap=cmap, norm=norm, linewidths=1.2)
    lc2d.set_array(vals)
    ax2d.add_collection(lc2d)
    xs = [p[0] for r in rows for p in (r["p0"], r["p1"])]
    ys = [p[1] for r in rows for p in (r["p0"], r["p1"])]
    ax2d.set_xlim(min(xs), max(xs))
    ax2d.set_ylim(min(ys), max(ys))
    ax2d.set_aspect("equal", adjustable="box")
    ax2d.set_xlabel("x (um)")
    ax2d.set_ylabel("y (um; + toward pia locally)")
    ax2d.set_title("x-y projection")

    cbar = fig.colorbar(lc2d, ax=[ax3d, ax2d], shrink=0.82)
    cbar.set_label("local irradiance at membrane (mW/mm2)")
    fig.suptitle(title)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def write_summary(summary, out: Path):
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--surface-irradiance-mw-mm2", type=float, default=1.0)
    parser.add_argument("--mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/light_distribution"))
    args = parser.parse_args()

    rule = load_pickle(args.cell_params)
    rows = collect_segments(rule)
    add_light(rows, args.surface_irradiance_mw_mm2, args.mu_eff_mm_inv)
    summary = summarize(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig_path = args.out_dir / f"PT5B_light_mu{args.mu_eff_mm_inv:g}_I{args.surface_irradiance_mw_mm2:g}.png"
    csv_path = args.out_dir / f"PT5B_light_summary_mu{args.mu_eff_mm_inv:g}_I{args.surface_irradiance_mw_mm2:g}.csv"
    plot(
        rows,
        fig_path,
        title=f"PT5B_full, 488 nm wide-field model, I0={args.surface_irradiance_mw_mm2:g} mW/mm2, mu_eff={args.mu_eff_mm_inv:g} mm^-1",
    )
    write_summary(summary, csv_path)
    for row in summary:
        print(
            f"{row['part']:16s} area={row['area_um2']:.1f} um2 "
            f"I_mean={row['area_weighted_irradiance_mw_mm2']:.4f} "
            f"I_range=[{row['min_irradiance_mw_mm2']:.4f}, {row['max_irradiance_mw_mm2']:.4f}] "
            f"power={row['intercepted_power_mw']:.6f} mW "
            f"fraction={row['fraction_of_non_axon_power']:.3f}"
        )
    print(fig_path)
    print(csv_path)


if __name__ == "__main__":
    main()
