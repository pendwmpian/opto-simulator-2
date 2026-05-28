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
from neuron import h

from calibrate_wang2007_gbar import (
    Wang2007Protocol,
    add_vitro_irradiance,
    current_metrics,
    half_max_x,
    segment_rows_vertical_slice,
    soma_axis_value,
)
from chr2_kinetics import ChR2FourState
from run_chr2_single_cell_sweep import CellFromNetPyNE


@dataclass(frozen=True)
class SEClampSuiteResult:
    label: str
    rs_mohm: float
    irradiance_scale: float
    gbar_mS_cm2: float
    peak_1s_9p2_nA: float
    time_to_peak_1s_9p2_ms: float
    inactivation_tau_1s_9p2_ms: float | None
    intensity_k_mw_mm2: float | None
    intensity_imax_nA: float
    duration_k_ms: float | None
    duration_imax_nA: float
    retained_area_um2: float
    mean_irradiance_mw_mm2: float


def soma_section(cell: CellFromNetPyNE):
    soma_name = "soma" if "soma" in cell.sections else next(name for name in cell.sections if name.startswith("soma"))
    return cell.sections[soma_name]


def prepare_rows(cell: CellFromNetPyNE, protocol: Wang2007Protocol, args):
    soma_z = soma_axis_value(cell, args.slice_axis)
    slice_center_um = soma_z + protocol.slice_thickness_um / 2.0 - args.soma_entry_depth_um
    illumination_entry_um = soma_z - args.soma_entry_depth_um
    rows, kept_fraction = segment_rows_vertical_slice(
        cell,
        args.slice_axis,
        slice_center_um,
        protocol.slice_thickness_um,
        args.truncate_radius_um,
    )
    add_vitro_irradiance(
        rows,
        protocol.irradiance_mw_mm2,
        args.vitro_mu_eff_mm_inv,
        "depth_decay",
        args.illumination_axis,
        illumination_entry_um,
    )
    for row in rows:
        row["base_irradiance_mw_mm2"] = row["irradiance_mw_mm2"]
    return rows, kept_fraction


def scaled_rows(rows, scale: float):
    return [
        {
            **row,
            "irradiance_mw_mm2": row.get("base_irradiance_mw_mm2", row["irradiance_mw_mm2"]) * scale,
        }
        for row in rows
    ]


def install_seclamp(cell: CellFromNetPyNE, rs_mohm: float, holding_mV: float, tstop_ms: float):
    clamp = h.SEClamp(soma_section(cell)(0.5))
    clamp.dur1 = tstop_ms + 1.0
    clamp.amp1 = holding_mV
    clamp.rs = rs_mohm
    return clamp


def install_chr2_current_sources(rows):
    clamps = []
    for row in rows:
        stim = h.IClamp(row["sec"](row["xloc"]))
        stim.delay = 0.0
        stim.dur = 1.0e9
        stim.amp = 0.0
        clamps.append(stim)
    return clamps


def simulate_seclamp_current(
    cell_params: Path,
    rows_template,
    gbar_mS_cm2: float,
    rs_mohm: float,
    protocol: Wang2007Protocol,
    dt_ms: float,
    tstop_ms: float,
    pulse_duration_ms: float,
    irradiance_scale: float = 1.0,
):
    h.load_file("stdrun.hoc")
    cell = CellFromNetPyNE(cell_params)

    # Re-map rows by section name and x location because NEURON sections are
    # recreated for every independent simulation.
    rows = []
    for template in rows_template:
        sec_name = template["sec"].name().split(".")[-1]
        sec = cell.sections[sec_name]
        rows.append({**template, "sec": sec})

    clamp = install_seclamp(cell, rs_mohm, protocol.holding_voltage_mV, tstop_ms)
    current_sources = install_chr2_current_sources(rows)
    chr2_models = [ChR2FourState() for _ in rows]

    h.dt = dt_ms
    h.tstop = tstop_ms
    h.finitialize(protocol.holding_voltage_mV)

    times = []
    currents = []
    while h.t <= tstop_ms:
        light_on = 0.0 <= h.t < pulse_duration_ms
        for row, source, model in zip(rows, current_sources, chr2_models):
            irr = row["irradiance_mw_mm2"] * irradiance_scale if light_on else 0.0
            model.step(dt_ms, irr)
            v_mV = float(row["sec"](row["xloc"]).v)
            inward_density_uA_cm2 = model.current_density_uA_cm2(gbar_mS_cm2, v_mV)
            source.amp = inward_density_uA_cm2 * row["area_cm2"] * 1000.0
        times.append(float(h.t))
        # Positive value means inward photocurrent that must be opposed by the clamp.
        currents.append(float(abs(clamp.i)))
        h.fadvance()
    return np.asarray(times), np.asarray(currents)


def calibrate_gbar_seclamp(cell_params, rows, rs_mohm, protocol, dt_ms, tstop_ms):
    lo = 1.0e-5
    hi = 1.0
    last = None
    for _ in range(16):
        mid = math.sqrt(lo * hi)
        t, current = simulate_seclamp_current(
            cell_params,
            rows,
            mid,
            rs_mohm,
            protocol,
            dt_ms,
            tstop_ms,
            pulse_duration_ms=1000.0,
        )
        peak = float(np.max(current))
        last = (mid, t, current, peak)
        if peak < protocol.target_peak_current_nA:
            lo = mid
        else:
            hi = mid
    if last is None:
        raise RuntimeError("SEClamp calibration did not run")
    gbar, t, current, _ = last
    return gbar, t, current


def run_suite_for_condition(label, cell_params, rows, rs_mohm, irradiance_scale, protocol, args):
    rows = scaled_rows(rows, irradiance_scale)
    gbar, t, current = calibrate_gbar_seclamp(cell_params, rows, rs_mohm, protocol, args.dt_ms, 1200.0)
    peak, _, time_to_peak, tau = current_metrics(t, current, protocol)

    intensity_values = [0.05, 0.1, 0.2, 0.4, 0.84, 1.0, 2.0, 4.6, 9.2]
    intensity_peaks = []
    for irr in intensity_values:
        _, current_i = simulate_seclamp_current(
            cell_params,
            rows,
            gbar,
            rs_mohm,
            protocol,
            args.dt_ms,
            180.0,
            pulse_duration_ms=100.0,
            irradiance_scale=irr / protocol.irradiance_mw_mm2,
        )
        intensity_peaks.append(float(np.max(current_i)))

    duration_values = [1, 2, 3, 4, 5, 8, 10, 20, 50, 100]
    duration_peaks = []
    for dur in duration_values:
        _, current_d = simulate_seclamp_current(
            cell_params,
            rows,
            gbar,
            rs_mohm,
            protocol,
            args.dt_ms,
            max(120.0, float(dur) + 60.0),
            pulse_duration_ms=float(dur),
        )
        duration_peaks.append(float(np.max(current_d)))

    total_area = sum(row["area_um2"] for row in rows)
    mean_irr = sum(row["irradiance_mw_mm2"] * row["area_um2"] for row in rows) / total_area
    result = SEClampSuiteResult(
        label=label,
        rs_mohm=rs_mohm,
        irradiance_scale=irradiance_scale,
        gbar_mS_cm2=float(gbar),
        peak_1s_9p2_nA=float(peak),
        time_to_peak_1s_9p2_ms=float(time_to_peak),
        inactivation_tau_1s_9p2_ms=tau,
        intensity_k_mw_mm2=half_max_x(intensity_values, intensity_peaks),
        intensity_imax_nA=float(max(intensity_peaks)),
        duration_k_ms=half_max_x([float(v) for v in duration_values], duration_peaks),
        duration_imax_nA=float(max(duration_peaks)),
        retained_area_um2=float(total_area),
        mean_irradiance_mw_mm2=float(mean_irr),
    )
    return result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks


def plot_condition(out, label, rs_mohm, result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks, protocol):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    axes[0].plot(t, current, color="#222222", lw=1.3)
    axes[0].axhline(protocol.target_peak_current_nA, color="#d62728", ls="--", lw=1)
    axes[0].set_xlim(-10, 250)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("SEClamp current magnitude (nA)")
    axes[0].set_title(f"1 s, 9.2 mW/mm2, Rs={rs_mohm:g}, scale={result.irradiance_scale:g}")

    axes[1].plot(intensity_values, intensity_peaks, marker="o")
    axes[1].axvline(0.84, color="#d62728", ls="--", lw=1)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("peak current (nA)")
    axes[1].set_title(f"intensity K={result.intensity_k_mw_mm2:.3g}")

    axes[2].plot(duration_values, duration_peaks, marker="o")
    axes[2].axvline(3.2, color="#d62728", ls="--", lw=1)
    axes[2].set_xscale("log")
    axes[2].set_xlabel("duration (ms)")
    axes[2].set_ylabel("peak current (nA)")
    axes[2].set_title(f"duration K={result.duration_k_ms:.3g}")
    fig.savefig(out / f"{label}_Rs{str(rs_mohm).replace('.', 'p')}_scale{str(result.irradiance_scale).replace('.', 'p')}_seclamp_suite.png", dpi=220)
    plt.close(fig)


def plot_scale_summary(results, out):
    by_label = {}
    for result in results:
        by_label.setdefault(result.label, []).append(result)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), constrained_layout=True)
    metrics = [
        ("time_to_peak_1s_9p2_ms", "time to peak (ms)", 12.0),
        ("intensity_k_mw_mm2", "intensity K (mW/mm2)", 0.84),
        ("duration_k_ms", "duration K (ms)", 3.2),
    ]
    for ax, (field, ylabel, target) in zip(axes, metrics):
        for label, values in by_label.items():
            values = sorted(values, key=lambda r: r.irradiance_scale)
            ax.plot([r.irradiance_scale for r in values], [getattr(r, field) for r in values], marker="o", label=label)
        ax.axhline(target, color="#d62728", ls="--", lw=1)
        ax.set_xlabel("effective irradiance scale")
        ax.set_ylabel(ylabel)
    axes[0].legend(frameon=False)
    fig.savefig(out / "seclamp_irradiance_scale_summary.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_seclamp_suite"))
    parser.add_argument("--pt-cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--it-cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/IT5B_full_cellParams.pkl"))
    parser.add_argument("--rs-mohm", default="5,7.5,10")
    parser.add_argument("--irradiance-scales", default="1.0")
    parser.add_argument("--soma-entry-depth-um", type=float, default=100.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--illumination-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    protocol = Wang2007Protocol()
    results = []
    for label, path in [("PT5B_full", args.pt_cell_params), ("IT5B_full", args.it_cell_params)]:
        h.load_file("stdrun.hoc")
        cell = CellFromNetPyNE(path)
        rows, _ = prepare_rows(cell, protocol, args)
        for rs in [float(v.strip()) for v in args.rs_mohm.split(",") if v.strip()]:
            for scale in [float(v.strip()) for v in args.irradiance_scales.split(",") if v.strip()]:
                result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks = run_suite_for_condition(
                    label, path, rows, rs, scale, protocol, args
                )
                results.append(result)
                plot_condition(args.out_dir, label, rs, result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks, protocol)
                print(result)

    with (args.out_dir / "seclamp_suite_results.json").open("w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)
    with (args.out_dir / "seclamp_suite_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in results)
    plot_scale_summary(results, args.out_dir)


if __name__ == "__main__":
    main()
