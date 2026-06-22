from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from neuron import h

from chr2_kinetics import ChR2FourState
from phase4_area_threshold import section_midpoint_xyz
from run_chr2_single_cell_sweep import CellFromNetPyNE, section_kind


@dataclass(frozen=True)
class Wang2007Protocol:
    preparation: str = "acute cortical parasagittal slice"
    line: str = "Thy1-ChR2-YFP line 18"
    cell_description: str = "ChR2-positive cortical layer V pyramidal neuron"
    slice_thickness_um: float = 300.0
    wavelength_nm_range: str = "465-495"
    large_field_area_mm2: float = 0.4
    irradiance_mw_mm2: float = 9.2
    pulse_duration_ms: float = 1000.0
    holding_voltage_mV: float = -70.0
    temperature_c: str = "21-24"
    blockers: str = "CNQX + APV + picrotoxin or GABAzine"
    target_peak_current_nA: float = 0.557
    target_peak_current_sd_nA: float = 0.153
    target_imax_nA: float = 0.642
    target_time_to_peak_ms: float = 12.0
    target_inactivation_tau_ms: float = 48.0


@dataclass(frozen=True)
class CalibrationResult:
    label: str
    cell_params: str
    gbar_mS_cm2: float
    peak_current_nA: float
    steady_current_nA: float
    time_to_peak_ms: float
    inactivation_tau_ms: float | None
    illuminated_area_um2: float
    mean_irradiance_mw_mm2: float
    total_membrane_area_um2: float
    slice_kept_fraction: float
    soma_entry_depth_um: float | None = None


@dataclass(frozen=True)
class ProtocolSuiteSummary:
    label: str
    gbar_mS_cm2: float
    peak_1s_9p2_nA: float
    time_to_peak_1s_9p2_ms: float
    inactivation_tau_1s_9p2_ms: float | None
    intensity_k_mw_mm2: float | None
    intensity_imax_nA: float
    duration_k_ms: float | None
    duration_imax_nA: float


def segment_rows(cell: CellFromNetPyNE, slice_depth_um: float, slice_thickness_um: float, truncate_radius_um: float | None):
    soma_name = "soma" if "soma" in cell.sections else next(name for name in cell.sections if name.startswith("soma"))
    soma = cell.sections[soma_name]
    soma_xyz = np.array([section_midpoint_xyz(soma, 0.5)])
    soma_y = float(soma_xyz[0, 1])
    rows = []
    total = 0
    kept = 0
    for sec_name, sec in cell.sections.items():
        kind = section_kind(sec_name)
        if kind == "axon":
            continue
        for seg in sec:
            total += 1
            x, y, z = section_midpoint_xyz(sec, seg.x)
            rel = np.array([x, y, z]) - soma_xyz[0]
            if truncate_radius_um is not None and float(np.linalg.norm(rel)) > truncate_radius_um:
                continue
            tissue_depth = slice_depth_um + (soma_y - y)
            if tissue_depth < 0 or tissue_depth > slice_thickness_um:
                continue
            area_um2 = float(h.area(seg.x, sec=sec))
            rows.append(
                {
                    "sec": sec,
                    "xloc": seg.x,
                    "kind": kind,
                    "x": x,
                    "y": y,
                    "z": z,
                    "tissue_depth_um": tissue_depth,
                    "area_um2": area_um2,
                    "area_cm2": area_um2 * 1.0e-8,
                }
            )
            kept += 1
    return rows, kept / total if total else 0.0


def segment_rows_vertical_slice(
    cell: CellFromNetPyNE,
    slice_axis: str,
    slice_center_um: float,
    slice_thickness_um: float,
    truncate_radius_um: float | None,
):
    soma_name = "soma" if "soma" in cell.sections else next(name for name in cell.sections if name.startswith("soma"))
    soma = cell.sections[soma_name]
    soma_xyz = np.array(section_midpoint_xyz(soma, 0.5), dtype=float)
    axis_index = {"x": 0, "y": 1, "z": 2}[slice_axis]
    half_thickness = slice_thickness_um / 2.0
    rows = []
    total = 0
    kept = 0
    for sec_name, sec in cell.sections.items():
        kind = section_kind(sec_name)
        if kind == "axon":
            continue
        for seg in sec:
            total += 1
            x, y, z = section_midpoint_xyz(sec, seg.x)
            xyz = np.array([x, y, z], dtype=float)
            rel = xyz - soma_xyz
            if truncate_radius_um is not None and float(np.linalg.norm(rel)) > truncate_radius_um:
                continue
            axis_value = float(xyz[axis_index])
            if abs(axis_value - slice_center_um) > half_thickness:
                continue
            area_um2 = float(h.area(seg.x, sec=sec))
            rows.append(
                {
                    "sec": sec,
                    "xloc": seg.x,
                    "kind": kind,
                    "x": x,
                    "y": y,
                    "z": z,
                    "slice_axis_value_um": axis_value,
                    "area_um2": area_um2,
                    "area_cm2": area_um2 * 1.0e-8,
                }
            )
            kept += 1
    return rows, kept / total if total else 0.0


def add_vitro_irradiance(
    rows,
    irradiance_mw_mm2: float,
    mu_eff_mm_inv: float,
    mode: str,
    illumination_axis: str = "z",
    illumination_entry_um: float | None = None,
):
    axis_key = {"x": "x", "y": "y", "z": "z"}[illumination_axis]
    if illumination_entry_um is None and rows:
        illumination_entry_um = min(float(row[axis_key]) for row in rows)
    for row in rows:
        if mode == "uniform":
            row["irradiance_mw_mm2"] = irradiance_mw_mm2
        elif mode == "depth_decay":
            depth_um = max(0.0, float(row[axis_key]) - float(illumination_entry_um))
            row["illumination_depth_um"] = depth_um
            depth_mm = depth_um / 1000.0
            row["irradiance_mw_mm2"] = irradiance_mw_mm2 * math.exp(-mu_eff_mm_inv * depth_mm)
        else:
            raise ValueError(f"Unknown vitro illumination mode: {mode}")


def simulate_current(
    rows,
    gbar_mS_cm2: float,
    protocol: Wang2007Protocol,
    dt_ms: float,
    tstop_ms: float,
    pulse_duration_ms: float | None = None,
    irradiance_scale: float = 1.0,
):
    models = [ChR2FourState() for _ in rows]
    pulse_duration = protocol.pulse_duration_ms if pulse_duration_ms is None else pulse_duration_ms
    t = np.arange(0.0, tstop_ms + dt_ms, dt_ms)
    current = np.zeros_like(t)
    for i, ti in enumerate(t):
        light_on = 0 <= ti < pulse_duration
        total_nA = 0.0
        for row, model in zip(rows, models):
            irr = row["irradiance_mw_mm2"] * irradiance_scale if light_on else 0.0
            model.step(dt_ms, irr)
            density_uA_cm2 = model.current_density_uA_cm2(gbar_mS_cm2, protocol.holding_voltage_mV)
            total_nA += density_uA_cm2 * row["area_cm2"] * 1000.0
        current[i] = total_nA
    return t, current


def current_metrics(t: np.ndarray, current: np.ndarray, protocol: Wang2007Protocol):
    peak_idx = int(np.argmax(current))
    peak = float(current[peak_idx])
    time_to_peak = float(t[peak_idx])
    steady_mask = (t >= protocol.pulse_duration_ms - 100.0) & (t < protocol.pulse_duration_ms)
    steady = float(np.mean(current[steady_mask])) if np.any(steady_mask) else float(current[-1])
    fit_mask = (t >= time_to_peak + 5.0) & (t <= min(protocol.pulse_duration_ms, time_to_peak + 250.0))
    tau = None
    if np.count_nonzero(fit_mask) > 10 and peak > steady:
        y = current[fit_mask] - steady
        ok = y > max(1e-9, 0.02 * (peak - steady))
        if np.count_nonzero(ok) > 10:
            x = t[fit_mask][ok] - time_to_peak
            log_y = np.log(y[ok])
            slope, _ = np.polyfit(x, log_y, 1)
            if slope < 0:
                tau = float(-1.0 / slope)
    return peak, steady, time_to_peak, tau


def calibrate_gbar(rows, protocol: Wang2007Protocol, dt_ms: float, tstop_ms: float, lo: float, hi: float, iterations: int):
    target = protocol.target_peak_current_nA
    last = None
    for _ in range(iterations):
        mid = math.sqrt(lo * hi)
        t, current = simulate_current(rows, mid, protocol, dt_ms, tstop_ms)
        peak, steady, time_to_peak, tau = current_metrics(t, current, protocol)
        last = (mid, t, current, peak, steady, time_to_peak, tau)
        if peak < target:
            lo = mid
        else:
            hi = mid
    return last


def calibrate_gbar_direct(rows, protocol: Wang2007Protocol, dt_ms: float, tstop_ms: float):
    unit_gbar = 1.0
    t, unit_current = simulate_current(rows, unit_gbar, protocol, dt_ms, tstop_ms)
    peak, _, _, _ = current_metrics(t, unit_current, protocol)
    if peak <= 0:
        raise ValueError("Unit-gbar current did not produce a positive peak current")
    gbar = protocol.target_peak_current_nA / peak
    current = unit_current * gbar
    peak, steady, time_to_peak, tau = current_metrics(t, current, protocol)
    return gbar, t, current, peak, steady, time_to_peak, tau


def half_max_x(xs: list[float], ys: list[float]) -> float | None:
    if not xs or not ys:
        return None
    ymin = min(ys)
    ymax = max(ys)
    target = ymin + 0.5 * (ymax - ymin)
    for i in range(1, len(xs)):
        y0, y1 = ys[i - 1], ys[i]
        if (y0 <= target <= y1) or (y1 <= target <= y0):
            if y1 == y0:
                return float(xs[i])
            frac = (target - y0) / (y1 - y0)
            return float(xs[i - 1] + frac * (xs[i] - xs[i - 1]))
    return None


def run_protocol_suite(
    rows,
    label: str,
    gbar_mS_cm2: float,
    protocol: Wang2007Protocol,
    dt_ms: float,
    out_dir: Path,
):
    t, current = simulate_current(rows, gbar_mS_cm2, protocol, dt_ms, 1200.0, pulse_duration_ms=1000.0)
    peak, _, time_to_peak, tau = current_metrics(t, current, protocol)

    intensity_values = [0.05, 0.1, 0.2, 0.4, 0.84, 1.0, 2.0, 4.6, 9.2]
    intensity_peaks = []
    for irr in intensity_values:
        t_i, current_i = simulate_current(
            rows,
            gbar_mS_cm2,
            protocol,
            dt_ms,
            180.0,
            pulse_duration_ms=100.0,
            irradiance_scale=irr / protocol.irradiance_mw_mm2,
        )
        intensity_peaks.append(float(np.max(current_i)))
    intensity_k = half_max_x(intensity_values, intensity_peaks)

    duration_values = [1, 2, 3, 4, 5, 8, 10, 20, 50, 100]
    duration_peaks = []
    for dur in duration_values:
        t_d, current_d = simulate_current(
            rows,
            gbar_mS_cm2,
            protocol,
            dt_ms,
            max(120.0, dur + 60.0),
            pulse_duration_ms=float(dur),
        )
        duration_peaks.append(float(np.max(current_d)))
    duration_k = half_max_x([float(v) for v in duration_values], duration_peaks)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    axes[0].plot(t, current, color="#222222", lw=1.3)
    axes[0].axhline(protocol.target_peak_current_nA, color="#d62728", ls="--", lw=1, label="557 pA target")
    axes[0].set_xlim(-10, 250)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("current (nA)")
    axes[0].set_title("1 s, 9.2 mW/mm2")
    axes[0].legend(frameon=False)

    axes[1].plot(intensity_values, intensity_peaks, marker="o", color="#1f77b4")
    axes[1].axvline(0.84, color="#d62728", ls="--", lw=1, label="Wang K=0.84")
    axes[1].axhline(protocol.target_imax_nA, color="#999999", ls=":", lw=1, label="Wang Imax=0.642")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("peak current (nA)")
    axes[1].set_title("100 ms intensity response")
    axes[1].legend(frameon=False)

    axes[2].plot(duration_values, duration_peaks, marker="o", color="#2ca02c")
    axes[2].axvline(3.2, color="#d62728", ls="--", lw=1, label="Wang K=3.2 ms")
    axes[2].set_xscale("log")
    axes[2].set_xlabel("pulse duration (ms)")
    axes[2].set_ylabel("peak current (nA)")
    axes[2].set_title("9.2 mW/mm2 duration response")
    axes[2].legend(frameon=False)

    fig.savefig(out_dir / f"{label}_wang_protocol_suite.png", dpi=220)
    plt.close(fig)

    with (out_dir / f"{label}_intensity_response.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["irradiance_mw_mm2", "peak_current_nA"])
        writer.writerows(zip(intensity_values, intensity_peaks))
    with (out_dir / f"{label}_duration_response.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pulse_duration_ms", "peak_current_nA"])
        writer.writerows(zip(duration_values, duration_peaks))

    return ProtocolSuiteSummary(
        label=label,
        gbar_mS_cm2=float(gbar_mS_cm2),
        peak_1s_9p2_nA=float(peak),
        time_to_peak_1s_9p2_ms=float(time_to_peak),
        inactivation_tau_1s_9p2_ms=tau,
        intensity_k_mw_mm2=intensity_k,
        intensity_imax_nA=float(max(intensity_peaks)),
        duration_k_ms=duration_k,
        duration_imax_nA=float(max(duration_peaks)),
    )


def soma_axis_value(cell: CellFromNetPyNE, axis: str) -> float:
    soma_name = "soma" if "soma" in cell.sections else next(name for name in cell.sections if name.startswith("soma"))
    soma = cell.sections[soma_name]
    xyz = section_midpoint_xyz(soma, 0.5)
    return float({"x": xyz[0], "y": xyz[1], "z": xyz[2]}[axis])


def plot_light_distribution(rows, label: str, out: Path, projection: str = "xy", view_from_light_entry: bool = False):
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    segs = []
    vals = []
    for row in rows:
        sec = row["sec"]
        xloc = float(row["xloc"])
        n3d = int(h.n3d(sec=sec))
        if n3d >= 2:
            target_arc = xloc * float(sec.L)
            half_len = float(sec.L) / max(1, int(sec.nseg)) / 2.0
            lo_arc = max(0.0, target_arc - half_len)
            hi_arc = min(float(sec.L), target_arc + half_len)
            p0 = point_at_arc(sec, lo_arc)
            p1 = point_at_arc(sec, hi_arc)
            if projection == "xy":
                segs.append([(p0[0], p0[1]), (p1[0], p1[1])])
            elif projection == "xz":
                segs.append([(p0[0], p0[2]), (p1[0], p1[2])])
            else:
                raise ValueError(projection)
        else:
            x, y, _ = section_midpoint_xyz(sec, xloc)
            dx = max(1.0, sec.L / max(1, sec.nseg)) / 2.0
            if projection == "xy":
                segs.append([(x - dx, y), (x + dx, y)])
            elif projection == "xz":
                _, _, z = section_midpoint_xyz(sec, xloc)
                segs.append([(x - dx, z), (x + dx, z)])
        vals.append(row["irradiance_mw_mm2"])
    lc = LineCollection(segs, cmap="viridis", linewidths=1.0)
    lc.set_array(np.asarray(vals))
    ax.add_collection(lc)
    xs = [p[0] for s in segs for p in s]
    ys = [p[1] for s in segs for p in s]
    ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(min(ys), max(ys))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    if projection == "xy":
        ax.set_ylabel("local morphology y (um)")
        if view_from_light_entry:
            ax.invert_xaxis()
            ax.set_title(f"{label}: x-y slice face viewed from light-entry side")
        else:
            ax.set_title(f"{label}: retained vertical slice x-y projection")
    else:
        ax.set_ylabel("slice thickness z (um)")
        if view_from_light_entry:
            # The model propagates light from low z toward high z. Put that
            # entry face at the top of the side-view figure for readability.
            ax.invert_yaxis()
            entry_row = min(rows, key=lambda row: row.get("illumination_depth_um", float("inf")))
            entry_z = float(entry_row["z"])
            soma_rows = [row for row in rows if row["kind"] == "soma"]
            soma_z = float(np.mean([row["z"] for row in soma_rows])) if soma_rows else float("nan")
            ax.axhline(entry_z, color="#d62728", ls="--", lw=1.2)
            ax.text(
                0.02,
                entry_z,
                " light-entry face",
                color="#d62728",
                fontsize=8,
                va="bottom",
                transform=ax.get_yaxis_transform(),
            )
            if math.isfinite(soma_z):
                ax.axhline(soma_z, color="#111111", ls=":", lw=1.0)
                ax.text(
                    0.98,
                    soma_z,
                    f"soma (~{soma_z - entry_z:.0f} um from entry) ",
                    color="#111111",
                    fontsize=8,
                    ha="right",
                    va="bottom",
                    transform=ax.get_yaxis_transform(),
                )
            ax.set_title(f"{label}: x-z side view (light enters from top)")
        else:
            ax.set_title(f"{label}: side view x-z projection")
    cbar = fig.colorbar(lc, ax=ax)
    cbar.set_label("irradiance (mW/mm2)")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_slice_light_schematic(rows, label: str, slice_thickness_um: float, illumination_axis: str, out: Path):
    fig, ax = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
    xs = [row["x"] for row in rows]
    zs = [row["z"] for row in rows]
    vals = [row["irradiance_mw_mm2"] for row in rows]
    sc = ax.scatter(xs, zs, c=vals, s=8, cmap="viridis")
    if zs:
        center_z = 0.5 * (min(zs) + max(zs))
    else:
        center_z = 0.0
    ax.axhspan(center_z - slice_thickness_um / 2, center_z + slice_thickness_um / 2, color="0.9", zorder=-10)
    entry_z = min(zs) if zs else -slice_thickness_um / 2
    ax.annotate(
        "light entry",
        xy=(min(xs) if xs else 0.0, entry_z),
        xytext=(min(xs) if xs else 0.0, entry_z - 80),
        arrowprops={"arrowstyle": "->", "color": "#d62728"},
        color="#d62728",
    )
    ax.set_xlabel("x (um)")
    ax.set_ylabel("z / slice thickness axis (um)")
    ax.set_title(f"{label}: vitro slice side view, light propagates along +{illumination_axis}")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("irradiance (mW/mm2)")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def point_at_arc(sec, target_arc: float) -> tuple[float, float, float]:
    n3d = int(h.n3d(sec=sec))
    if n3d == 0:
        return 0.0, 0.0, 0.0
    if n3d == 1:
        return float(h.x3d(0, sec=sec)), float(h.y3d(0, sec=sec)), float(h.z3d(0, sec=sec))
    prev_arc = float(h.arc3d(0, sec=sec))
    prev = np.array([h.x3d(0, sec=sec), h.y3d(0, sec=sec), h.z3d(0, sec=sec)], dtype=float)
    for i in range(1, n3d):
        arc = float(h.arc3d(i, sec=sec))
        cur = np.array([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)], dtype=float)
        if arc >= target_arc:
            frac = min(1.0, max(0.0, (target_arc - prev_arc) / max(1.0e-9, arc - prev_arc)))
            p = prev + frac * (cur - prev)
            return float(p[0]), float(p[1]), float(p[2])
        prev_arc = arc
        prev = cur
    return float(prev[0]), float(prev[1]), float(prev[2])


def plot_current(t: np.ndarray, current: np.ndarray, label: str, result: CalibrationResult, protocol: Wang2007Protocol, out: Path):
    fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
    ax.plot(t, current, color="#222222", lw=1.4)
    ax.axhline(protocol.target_peak_current_nA, color="#d62728", ls="--", lw=1, label="Wang mean peak")
    ax.fill_between(
        t,
        protocol.target_peak_current_nA - protocol.target_peak_current_sd_nA,
        protocol.target_peak_current_nA + protocol.target_peak_current_sd_nA,
        color="#d62728",
        alpha=0.12,
        label="Wang mean +/- SD",
    )
    ax.set_xlim(-10, min(1100, t[-1]))
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("total ChR2 inward current (nA)")
    ax.set_title(
        f"{label}: gbar={result.gbar_mS_cm2:.5g} mS/cm2, "
        f"peak={result.peak_current_nA:.3f} nA, tau={result.inactivation_tau_ms}"
    )
    ax.legend(frameon=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_soma_depth_sweep(results: list[CalibrationResult], out: Path):
    by_label = {}
    for result in results:
        by_label.setdefault(result.label, []).append(result)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    metrics = [
        ("gbar_mS_cm2", "gbar (mS/cm2)"),
        ("mean_irradiance_mw_mm2", "area-weighted irradiance (mW/mm2)"),
        ("time_to_peak_ms", "time to peak (ms)"),
    ]
    for ax, (field, ylabel) in zip(axes, metrics):
        for label, values in by_label.items():
            values = sorted(values, key=lambda r: r.soma_entry_depth_um or 0.0)
            xs = [r.soma_entry_depth_um for r in values]
            ys = [getattr(r, field) for r in values]
            ax.plot(xs, ys, marker="o", label=label)
        ax.axvspan(50, 100, color="0.9", zorder=-10, label="common target range" if field == "gbar_mS_cm2" else None)
        ax.set_xlabel("soma depth from light-entry surface (um)")
        ax.set_ylabel(ylabel)
    axes[0].legend(frameon=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def run_cell(label: str, cell_params: Path, args, protocol: Wang2007Protocol, soma_entry_depth_um: float | None = None):
    h.load_file("stdrun.hoc")
    cell = CellFromNetPyNE(cell_params)
    slice_center_um = args.slice_center_um
    illumination_entry_um = args.illumination_entry_um
    if soma_entry_depth_um is not None:
        soma_axis = soma_axis_value(cell, args.slice_axis)
        slice_center_um = soma_axis + protocol.slice_thickness_um / 2.0 - soma_entry_depth_um
        illumination_entry_um = soma_axis - soma_entry_depth_um
    if args.slice_model == "vertical":
        rows, kept_fraction = segment_rows_vertical_slice(
            cell,
            args.slice_axis,
            slice_center_um,
            protocol.slice_thickness_um,
            args.truncate_radius_um,
        )
    else:
        rows, kept_fraction = segment_rows(cell, args.slice_depth_um, protocol.slice_thickness_um, args.truncate_radius_um)
    add_vitro_irradiance(
        rows,
        protocol.irradiance_mw_mm2,
        args.vitro_mu_eff_mm_inv,
        args.vitro_illumination_mode,
        args.illumination_axis,
        illumination_entry_um,
    )
    if args.calibration_method == "direct":
        gbar, t, current, peak, steady, time_to_peak, tau = calibrate_gbar_direct(
            rows,
            protocol,
            args.dt_ms,
            args.tstop_ms,
        )
    else:
        gbar, t, current, peak, steady, time_to_peak, tau = calibrate_gbar(
            rows,
            protocol,
            args.dt_ms,
            args.tstop_ms,
            args.gbar_lo,
            args.gbar_hi,
            args.binary_iterations,
        )
    total_area = sum(row["area_um2"] for row in rows)
    mean_irr = sum(row["irradiance_mw_mm2"] * row["area_um2"] for row in rows) / total_area
    result = CalibrationResult(
        label=label,
        cell_params=str(cell_params),
        gbar_mS_cm2=float(gbar),
        peak_current_nA=float(peak),
        steady_current_nA=float(steady),
        time_to_peak_ms=float(time_to_peak),
        inactivation_tau_ms=tau,
        illuminated_area_um2=float(total_area),
        mean_irradiance_mw_mm2=float(mean_irr),
        total_membrane_area_um2=float(total_area),
        slice_kept_fraction=float(kept_fraction),
        soma_entry_depth_um=soma_entry_depth_um,
    )
    return result, rows, t, current


def parse_depths(raw: str) -> list[float]:
    return [float(part.strip()) for part in raw.split(",") if part.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_gbar_calibration"))
    parser.add_argument("--pt-cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--it-cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/IT5B_full_cellParams.pkl"))
    parser.add_argument("--slice-depth-um", type=float, default=150.0)
    parser.add_argument("--slice-model", choices=["vertical", "depth_window"], default="vertical")
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--slice-center-um", type=float, default=0.0)
    parser.add_argument("--illumination-axis", choices=["x", "z"], default="z")
    parser.add_argument("--illumination-entry-um", type=float, default=None)
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--vitro-illumination-mode", choices=["uniform", "depth_decay"], default="depth_decay")
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--dt-ms", type=float, default=0.05)
    parser.add_argument("--tstop-ms", type=float, default=1200.0)
    parser.add_argument("--gbar-lo", type=float, default=1e-5)
    parser.add_argument("--gbar-hi", type=float, default=100.0)
    parser.add_argument("--binary-iterations", type=int, default=20)
    parser.add_argument("--calibration-method", choices=["direct", "binary"], default="direct")
    parser.add_argument("--soma-entry-depths-um", default="50,75,100,150,200,250")
    args = parser.parse_args()

    protocol = Wang2007Protocol()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    protocol_suite_results = []
    for label, path in [("PT5B_full", args.pt_cell_params), ("IT5B_full", args.it_cell_params)]:
        result, rows, t, current = run_cell(label, path, args, protocol)
        results.append(result)
        plot_light_distribution(
            rows,
            label,
            args.out_dir / f"{label}_vitro_irradiance_xy.png",
            projection="xy",
            view_from_light_entry=True,
        )
        plot_light_distribution(
            rows,
            label,
            args.out_dir / f"{label}_vitro_irradiance_xz.png",
            projection="xz",
            view_from_light_entry=True,
        )
        plot_slice_light_schematic(rows, label, protocol.slice_thickness_um, args.illumination_axis, args.out_dir / f"{label}_slice_light_sideview.png")
        plot_current(t, current, label, result, protocol, args.out_dir / f"{label}_calibrated_current.png")
        protocol_suite_results.append(run_protocol_suite(rows, label, result.gbar_mS_cm2, protocol, args.dt_ms, args.out_dir))
        print(result)

    with (args.out_dir / "wang2007_protocol.json").open("w") as f:
        json.dump(asdict(protocol) | {"simulation_args": vars(args) | {"out_dir": str(args.out_dir), "pt_cell_params": str(args.pt_cell_params), "it_cell_params": str(args.it_cell_params)}}, f, indent=2)
    with (args.out_dir / "calibration_results.json").open("w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)
    with (args.out_dir / "calibration_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in results)
    with (args.out_dir / "wang_protocol_suite_results.json").open("w") as f:
        json.dump([asdict(r) for r in protocol_suite_results], f, indent=2)
    with (args.out_dir / "wang_protocol_suite_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(protocol_suite_results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in protocol_suite_results)

    sweep_results = []
    for depth_um in parse_depths(args.soma_entry_depths_um):
        for label, path in [("PT5B_full", args.pt_cell_params), ("IT5B_full", args.it_cell_params)]:
            result, _, _, _ = run_cell(label, path, args, protocol, soma_entry_depth_um=depth_um)
            sweep_results.append(result)
    with (args.out_dir / "soma_entry_depth_sweep.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(sweep_results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in sweep_results)
    with (args.out_dir / "soma_entry_depth_sweep.json").open("w") as f:
        json.dump([asdict(r) for r in sweep_results], f, indent=2)
    plot_soma_depth_sweep(sweep_results, args.out_dir / "soma_entry_depth_sweep.png")


if __name__ == "__main__":
    main()
