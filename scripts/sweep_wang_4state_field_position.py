from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from neuron import h

from calibrate_wang2007_gbar import Wang2007Protocol, point_at_arc
from rerun_wang_4state_no_truncate import (
    area_summary,
    build_no_truncate_rows,
    field_axes,
    run_trace,
)
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    calibrate_gbar,
    fit_hill_k,
    load_mechanisms,
    williams_q10_scales,
)


def parse_offsets(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def weighted_center(rows: list[dict], axis_a: str, axis_b: str, mode: str) -> tuple[float, float]:
    if mode == "soma":
        selected = [row for row in rows if row["kind"] == "soma"]
    elif mode == "apical":
        selected = [row for row in rows if row["kind"] == "apic"]
    elif mode == "all":
        selected = rows
    else:
        raise ValueError(f"Unknown center mode: {mode}")
    if not selected:
        selected = rows
    weights = np.asarray([row["area_um2"] for row in selected], dtype=float)
    center_a = float(np.average([row[axis_a] for row in selected], weights=weights))
    center_b = float(np.average([row[axis_b] for row in selected], weights=weights))
    return center_a, center_b


def apply_circular_field(
    rows: list[dict],
    center_a: float,
    center_b: float,
    radius_um: float,
    axis_a: str,
    axis_b: str,
) -> list[dict]:
    copied = [{**row, "localization_weight": 1.0, "field_weight": 1.0} for row in rows]
    for row in copied:
        dist = math.hypot(float(row[axis_a]) - center_a, float(row[axis_b]) - center_b)
        if dist > radius_um:
            row["irradiance_mw_mm2"] = 0.0
            row["field_weight"] = 0.0
    return copied


def make_positions(rows: list[dict], args) -> list[dict]:
    axis_a, axis_b = field_axes(args.illumination_axis)
    base_a, base_b = weighted_center(rows, axis_a, axis_b, args.center_mode)
    offsets = parse_offsets(args.offsets_um)
    if args.position_pattern == "cross":
        nonzero = [offset for offset in offsets if abs(offset) > 1.0e-9]
        step = nonzero[0] if nonzero else 200.0
        pairs = [(0.0, 0.0), (-step, 0.0), (step, 0.0), (0.0, -step), (0.0, step)]
    else:
        pairs = [(da, db) for da in offsets for db in offsets]
    positions = []
    for da, db in pairs:
        positions.append(
            {
                "position": f"{args.center_mode}_da{da:+.0f}_db{db:+.0f}",
                "center_a_um": base_a + da,
                "center_b_um": base_b + db,
                "offset_a_um": da,
                "offset_b_um": db,
            }
        )
    return positions


def simulate_position(args, protocol, base_rows: list[dict], position: dict, radius_um: float):
    axis_a, axis_b = field_axes(args.illumination_axis)
    rows = apply_circular_field(base_rows, position["center_a_um"], position["center_b_um"], radius_um, axis_a, axis_b)
    q10_scales = williams_q10_scales(args.reference_c, args.target_c)
    gbar, _, _ = calibrate_gbar(
        args.cell_params,
        rows,
        args.irradiance_scale,
        args.rs_mohm,
        protocol,
        args.dt_ms,
        args.binary_iterations,
        q10_scales,
        True,
        protocol.target_imax_nA,
        args.saturating_irradiance_mw_mm2,
        args.saturating_duration_ms,
        args.use_original_biophysics,
    )

    trace_1s = run_trace(args.cell_params, rows, gbar, protocol.irradiance_mw_mm2, 1000.0, 1200.0, args, protocol, q10_scales)
    from calibrate_wang2007_gbar import current_metrics

    peak, steady, ttp, tau = current_metrics(trace_1s["t"], trace_1s["filtered_i"], protocol)

    intensity_rows = []
    intensity_peaks = []
    for irr in WANG_INTENSITY_VALUES_MW_MM2:
        trace = run_trace(args.cell_params, rows, gbar, irr, 100.0, 180.0, args, protocol, q10_scales)
        p_i, s_i, ttp_i, tau_i = current_metrics(trace["t"], trace["filtered_i"], protocol)
        intensity_peaks.append(p_i)
        intensity_rows.append(
            {
                "position": position["position"],
                "irradiance_mw_mm2": irr,
                "filtered_peak_nA": p_i,
                "filtered_steady_nA": s_i,
                "filtered_ttp_ms": ttp_i,
                "filtered_tau_ms": tau_i,
            }
        )
    intensity_k, intensity_imax, intensity_hill = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, intensity_peaks)

    duration_values = [1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    duration_rows = []
    duration_peaks = []
    for dur in duration_values:
        trace = run_trace(args.cell_params, rows, gbar, protocol.irradiance_mw_mm2, float(dur), max(140.0, float(dur) + 80.0), args, protocol, q10_scales)
        p_d, s_d, ttp_d, tau_d = current_metrics(trace["t"], trace["filtered_i"], protocol)
        duration_peaks.append(p_d)
        duration_rows.append(
            {
                "position": position["position"],
                "duration_ms": float(dur),
                "filtered_peak_nA": p_d,
                "filtered_steady_nA": s_d,
                "filtered_ttp_ms": ttp_d,
                "filtered_tau_ms": tau_d,
            }
        )
    duration_k, duration_imax, duration_hill = fit_hill_k([float(v) for v in duration_values], duration_peaks)

    summary = {
        **position,
        "field_radius_um": radius_um,
        "gbar_mS_cm2": float(gbar),
        "filtered_peak_9p2_1s_nA": peak,
        "filtered_steady_9p2_1s_nA": steady,
        "filtered_steady_peak_ratio": steady / peak if peak else float("nan"),
        "filtered_ttp_9p2_1s_ms": ttp,
        "filtered_tau_9p2_1s_ms": tau,
        "filtered_intensity_k_mw_mm2": intensity_k,
        "filtered_intensity_imax_nA": intensity_imax,
        "filtered_intensity_hill_n": intensity_hill,
        "filtered_duration_k_ms": duration_k,
        "filtered_duration_imax_nA": duration_imax,
        "filtered_duration_hill_n": duration_hill,
        **area_summary(rows),
    }
    return rows, summary, intensity_rows, duration_rows, trace_1s, WANG_INTENSITY_VALUES_MW_MM2, intensity_peaks, duration_values, duration_peaks


def morphology_segments(rows: list[dict], axis_a: str, axis_b: str):
    segs = []
    colors = []
    for row in rows:
        sec = row["sec"]
        xloc = float(row["xloc"])
        half_len = float(sec.L) / max(1, int(sec.nseg)) / 2.0
        lo_arc = max(0.0, xloc * float(sec.L) - half_len)
        hi_arc = min(float(sec.L), xloc * float(sec.L) + half_len)
        p0 = point_at_arc(sec, lo_arc)
        p1 = point_at_arc(sec, hi_arc)
        idx = {"x": 0, "y": 1, "z": 2}
        segs.append([(p0[idx[axis_a]], p0[idx[axis_b]]), (p1[idx[axis_a]], p1[idx[axis_b]])])
        colors.append(row["irradiance_mw_mm2"])
    return segs, np.asarray(colors, dtype=float)


def plot_composite(
    args,
    protocol,
    rows: list[dict],
    summary: dict,
    trace_1s: dict,
    intensity_values,
    intensity_peaks,
    duration_values,
    duration_peaks,
    out: Path,
):
    axis_a, axis_b = field_axes(args.illumination_axis)
    fig = plt.figure(figsize=(15.5, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.0, 1.0])
    ax_map = fig.add_subplot(gs[:, 0])
    ax_a = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_txt = fig.add_subplot(gs[1, 2])

    segs, colors = morphology_segments(rows, axis_a, axis_b)
    lc = LineCollection(segs, cmap="viridis", linewidths=1.1)
    lc.set_array(colors)
    ax_map.add_collection(lc)
    xs = [p[0] for seg in segs for p in seg]
    ys = [p[1] for seg in segs for p in seg]
    circle = plt.Circle(
        (summary["center_a_um"], summary["center_b_um"]),
        summary["field_radius_um"],
        facecolor="none",
        edgecolor="#1f77ff",
        lw=2.2,
        alpha=0.95,
    )
    ax_map.add_patch(circle)
    ax_map.scatter([summary["center_a_um"]], [summary["center_b_um"]], s=24, color="#1f77ff", zorder=5)
    margin = summary["field_radius_um"] * 0.25
    ax_map.set_xlim(min(min(xs), summary["center_a_um"] - summary["field_radius_um"]) - margin, max(max(xs), summary["center_a_um"] + summary["field_radius_um"]) + margin)
    ax_map.set_ylim(min(min(ys), summary["center_b_um"] - summary["field_radius_um"]) - margin, max(max(ys), summary["center_b_um"] + summary["field_radius_um"]) + margin)
    ax_map.set_aspect("equal", adjustable="box")
    ax_map.set_xlabel(f"{axis_a} (um)")
    ax_map.set_ylabel(f"{axis_b} (um)")
    ax_map.set_title("0.4 mm2 circular field on slice face")
    cbar = fig.colorbar(lc, ax=ax_map, fraction=0.046, pad=0.04)
    cbar.set_label("local irradiance (mW/mm2)")

    ax_a.plot(trace_1s["t"], -trace_1s["filtered_i"], color="#222222", lw=1.2)
    ax_a.axhline(-protocol.target_peak_current_nA, color="#d62728", ls="--", lw=1, label="Wang peak")
    ax_a.set_xlim(-10, 250)
    ax_a.set_xlabel("time (ms)")
    ax_a.set_ylabel("current (nA)")
    ax_a.set_title("Fig. 2A-like: 9.2 mW/mm2, 1 s")
    ax_a.legend(fontsize=8, frameon=False)

    ax_c.plot(intensity_values, intensity_peaks, marker="o", color="#1f77b4")
    ax_c.axvline(0.84, color="#d62728", ls="--", lw=1, label="Wang K")
    ax_c.axhline(protocol.target_imax_nA, color="#777777", ls=":", lw=1, label="Wang Imax")
    ax_c.set_xscale("log")
    ax_c.set_xlabel("irradiance (mW/mm2)")
    ax_c.set_ylabel("peak current (nA)")
    ax_c.set_title("Fig. 2C-like: intensity response")
    ax_c.legend(fontsize=8, frameon=False)

    ax_d.plot(duration_values, duration_peaks, marker="o", color="#2ca02c")
    ax_d.axvline(3.2, color="#d62728", ls="--", lw=1, label="Wang duration K")
    ax_d.set_xscale("log")
    ax_d.set_xlabel("pulse duration (ms)")
    ax_d.set_ylabel("peak current (nA)")
    ax_d.set_title("Fig. 2D/F-like: duration response")
    ax_d.legend(fontsize=8, frameon=False)

    ax_txt.axis("off")
    text = (
        f"{summary['position']}\n"
        f"center ({axis_a}, {axis_b}) = ({summary['center_a_um']:.1f}, {summary['center_b_um']:.1f}) um\n"
        f"radius = {summary['field_radius_um']:.1f} um; area = {protocol.large_field_area_mm2:.2f} mm2\n\n"
        f"gbar = {summary['gbar_mS_cm2']:.5g} mS/cm2\n"
        f"peak 9.2/1s = {summary['filtered_peak_9p2_1s_nA']:.3f} nA\n"
        f"steady/peak = {summary['filtered_steady_peak_ratio']:.3f}\n"
        f"TTP = {summary['filtered_ttp_9p2_1s_ms']:.2f} ms\n"
        f"decay tau = {summary['filtered_tau_9p2_1s_ms']:.1f} ms\n"
        f"K intensity = {summary['filtered_intensity_k_mw_mm2']:.3g} mW/mm2\n"
        f"K duration = {summary['filtered_duration_k_ms']:.3g} ms\n\n"
        f"retained area = {summary['retained_area_um2']:.0f} um2\n"
        f"illuminated area = {summary['illuminated_area_um2']:.0f} um2\n"
        f"illuminated apic frac = {summary['illuminated_apic_fraction']:.3f}"
    )
    ax_txt.text(0.0, 1.0, text, va="top", ha="left", family="monospace", fontsize=10)

    fig.suptitle("Wang 2007 4-state 22C field-position sweep, no morphology truncation", fontsize=13)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_sweep_summary(summaries: list[dict], out: Path):
    if not summaries:
        return
    x = [row["center_a_um"] for row in summaries]
    y = [row["center_b_um"] for row in summaries]
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 3.8), constrained_layout=True)
    metrics = [
        ("gbar_mS_cm2", "gbar (mS/cm2)"),
        ("filtered_ttp_9p2_1s_ms", "TTP (ms)"),
        ("filtered_intensity_k_mw_mm2", "intensity K"),
        ("filtered_duration_k_ms", "duration K"),
    ]
    for ax, (key, title) in zip(axes, metrics):
        vals = [row[key] for row in summaries]
        sc = ax.scatter(x, y, c=vals, s=110, cmap="viridis")
        for row in summaries:
            ax.text(row["center_a_um"], row["center_b_um"], row["position"].split("_", 1)[1], fontsize=6, ha="center", va="center", color="white")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("field center axis-a (um)")
        ax.set_ylabel("field center axis-b (um)")
        ax.set_title(title)
        fig.colorbar(sc, ax=ax)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang_4state_field_position_sweep"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--center-mode", choices=["apical", "soma", "all"], default="apical")
    parser.add_argument("--offsets-um", default="-200,0,200")
    parser.add_argument("--position-pattern", choices=["grid", "cross"], default="grid")
    parser.add_argument("--rs-mohm", type=float, default=10.0)
    parser.add_argument("--electrode-rs-mohm", type=float, default=10.0)
    parser.add_argument("--pipette-capacitance-pf", type=float, default=100.0)
    parser.add_argument("--amplifier-filter-tau-ms", type=float, default=0.25)
    parser.add_argument("--amplifier-filter-order", type=int, default=4)
    parser.add_argument("--irradiance-scale", type=float, default=1.0)
    parser.add_argument("--soma-entry-depth-um", type=float, default=100.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--illumination-axis", choices=["x", "z"], default="z")
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=1.3)
    parser.add_argument("--reference-c", type=float, default=37.0)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--dt-ms", type=float, default=0.05)
    parser.add_argument("--binary-iterations", type=int, default=10)
    parser.add_argument("--saturating-irradiance-mw-mm2", type=float, default=100.0)
    parser.add_argument("--saturating-duration-ms", type=float, default=100.0)
    parser.add_argument("--use-original-biophysics", action="store_true", default=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    protocol = Wang2007Protocol()
    h.load_file("stdrun.hoc")
    cell = CellFromNetPyNE(args.cell_params)
    base_rows = build_no_truncate_rows(cell, protocol, args)
    radius_um = math.sqrt(protocol.large_field_area_mm2 * 1.0e6 / math.pi)
    positions = make_positions(base_rows, args)

    summaries = []
    all_intensity_rows = []
    all_duration_rows = []
    for position in positions:
        print(f"running {position['position']}")
        rows, summary, intensity_rows, duration_rows, trace_1s, intensity_values, intensity_peaks, duration_values, duration_peaks = simulate_position(
            args, protocol, base_rows, position, radius_um
        )
        summaries.append(summary)
        all_intensity_rows.extend(intensity_rows)
        all_duration_rows.extend(duration_rows)
        plot_composite(
            args,
            protocol,
            rows,
            summary,
            trace_1s,
            intensity_values,
            intensity_peaks,
            duration_values,
            duration_peaks,
            args.out_dir / f"{position['position']}_composite.png",
        )

    with (args.out_dir / "field_position_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(summaries)
    with (args.out_dir / "field_position_intensity.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_intensity_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_intensity_rows)
    with (args.out_dir / "field_position_duration.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_duration_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_duration_rows)
    plot_sweep_summary(summaries, args.out_dir / "field_position_sweep_summary.png")


if __name__ == "__main__":
    main()
