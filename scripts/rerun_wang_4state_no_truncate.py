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
from neuron import h

from calibrate_wang2007_gbar import (
    Wang2007Protocol,
    add_vitro_irradiance,
    current_metrics,
    segment_rows_vertical_slice,
    soma_axis_value,
)
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_electrode_amplifier_model import apply_electrode_amplifier_readout
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    calibrate_gbar,
    fit_hill_k,
    load_mechanisms,
    simulate_mod_current,
    williams_q10_scales,
)


def field_axes(illumination_axis: str) -> tuple[str, str]:
    axes = ["x", "y", "z"]
    axes.remove(illumination_axis)
    return axes[0], axes[1]


def build_no_truncate_rows(cell: CellFromNetPyNE, protocol: Wang2007Protocol, args):
    soma_axis = soma_axis_value(cell, args.slice_axis)
    slice_center_um = soma_axis + protocol.slice_thickness_um / 2.0 - args.soma_entry_depth_um
    illumination_entry_um = soma_axis - args.soma_entry_depth_um
    rows, _ = segment_rows_vertical_slice(
        cell,
        args.slice_axis,
        slice_center_um,
        protocol.slice_thickness_um,
        None,
    )
    add_vitro_irradiance(
        rows,
        protocol.irradiance_mw_mm2,
        args.vitro_mu_eff_mm_inv,
        "depth_decay",
        args.illumination_axis,
        illumination_entry_um,
    )
    return rows


def apply_illumination_mode(rows: list[dict], mode: str, field_area_mm2: float, illumination_axis: str):
    copied = [{**row, "localization_weight": 1.0, "field_weight": 1.0} for row in rows]
    if mode == "full":
        return copied, {
            "field_center_a_um": float("nan"),
            "field_center_b_um": float("nan"),
            "field_radius_um": float("nan"),
        }
    if mode != "apical_0p4mm2":
        raise ValueError(f"Unknown illumination mode: {mode}")

    axis_a, axis_b = field_axes(illumination_axis)
    apical = [row for row in copied if row["kind"] == "apic"]
    if not apical:
        apical = copied
    area = np.asarray([row["area_um2"] for row in apical], dtype=float)
    center_a = float(np.average([row[axis_a] for row in apical], weights=area))
    center_b = float(np.average([row[axis_b] for row in apical], weights=area))
    radius_um = math.sqrt(field_area_mm2 * 1.0e6 / math.pi)

    for row in copied:
        dist = math.hypot(float(row[axis_a]) - center_a, float(row[axis_b]) - center_b)
        if dist > radius_um:
            row["irradiance_mw_mm2"] = 0.0
            row["field_weight"] = 0.0
    return copied, {
        "field_center_a_um": center_a,
        "field_center_b_um": center_b,
        "field_radius_um": radius_um,
    }


def area_summary(rows: list[dict]) -> dict:
    retained_area = float(sum(row["area_um2"] for row in rows))
    illuminated_area = float(sum(row["area_um2"] for row in rows if row["irradiance_mw_mm2"] > 0.0))
    mean_irr = (
        float(sum(row["irradiance_mw_mm2"] * row["area_um2"] for row in rows) / retained_area)
        if retained_area
        else float("nan")
    )
    apic_area = float(sum(row["area_um2"] for row in rows if row["kind"] == "apic"))
    illum_apic_area = float(sum(row["area_um2"] for row in rows if row["kind"] == "apic" and row["irradiance_mw_mm2"] > 0.0))
    return {
        "retained_area_um2": retained_area,
        "illuminated_area_um2": illuminated_area,
        "illuminated_area_fraction": illuminated_area / retained_area if retained_area else float("nan"),
        "apic_area_um2": apic_area,
        "illuminated_apic_area_um2": illum_apic_area,
        "illuminated_apic_fraction": illum_apic_area / apic_area if apic_area else float("nan"),
        "mean_irradiance_mw_mm2": mean_irr,
    }


def run_trace(cell_params, rows, gbar, irradiance, pulse_duration_ms, tstop_ms, args, protocol, q10_scales):
    t, clamp_i, soma_v = simulate_mod_current(
        cell_params,
        rows,
        gbar,
        args.irradiance_scale,
        args.rs_mohm,
        protocol,
        args.dt_ms,
        tstop_ms,
        pulse_duration_ms,
        protocol_irradiance_scale=irradiance / protocol.irradiance_mw_mm2,
        q10_scales=q10_scales,
        use_photon_flux=True,
        use_original_biophysics=args.use_original_biophysics,
        subtract_no_light=args.use_original_biophysics,
    )
    return {"t": t, "clamp_i": clamp_i, "filtered_i": apply_electrode_amplifier_readout(t, clamp_i, args), "soma_v": soma_v}


def run_mode(args, protocol, base_rows, mode: str):
    rows, field = apply_illumination_mode(base_rows, mode, protocol.large_field_area_mm2, args.illumination_axis)
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
    raw_peak, raw_steady, raw_ttp, raw_tau = current_metrics(trace_1s["t"], trace_1s["clamp_i"], protocol)
    filt_peak, filt_steady, filt_ttp, filt_tau = current_metrics(trace_1s["t"], trace_1s["filtered_i"], protocol)

    intensity_rows = []
    filtered_peaks = []
    raw_peaks = []
    for irr in WANG_INTENSITY_VALUES_MW_MM2:
        trace = run_trace(args.cell_params, rows, gbar, irr, 100.0, 180.0, args, protocol, q10_scales)
        raw_peak_i, raw_steady_i, raw_ttp_i, raw_tau_i = current_metrics(trace["t"], trace["clamp_i"], protocol)
        filt_peak_i, filt_steady_i, filt_ttp_i, filt_tau_i = current_metrics(trace["t"], trace["filtered_i"], protocol)
        raw_peaks.append(raw_peak_i)
        filtered_peaks.append(filt_peak_i)
        intensity_rows.append(
            {
                "mode": mode,
                "irradiance_mw_mm2": irr,
                "raw_peak_nA": raw_peak_i,
                "raw_steady_nA": raw_steady_i,
                "raw_ttp_ms": raw_ttp_i,
                "raw_tau_ms": raw_tau_i,
                "filtered_peak_nA": filt_peak_i,
                "filtered_steady_nA": filt_steady_i,
                "filtered_ttp_ms": filt_ttp_i,
                "filtered_tau_ms": filt_tau_i,
            }
        )

    raw_k, raw_imax, raw_hill = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, raw_peaks)
    filt_k, filt_imax, filt_hill = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, filtered_peaks)
    summary = {
        "mode": mode,
        "gbar_mS_cm2": gbar,
        "raw_peak_9p2_1s_nA": raw_peak,
        "raw_steady_9p2_1s_nA": raw_steady,
        "raw_steady_peak_ratio": raw_steady / raw_peak if raw_peak else float("nan"),
        "raw_ttp_9p2_1s_ms": raw_ttp,
        "raw_tau_9p2_1s_ms": raw_tau,
        "filtered_peak_9p2_1s_nA": filt_peak,
        "filtered_steady_9p2_1s_nA": filt_steady,
        "filtered_steady_peak_ratio": filt_steady / filt_peak if filt_peak else float("nan"),
        "filtered_ttp_9p2_1s_ms": filt_ttp,
        "filtered_tau_9p2_1s_ms": filt_tau,
        "raw_intensity_k_mw_mm2": raw_k,
        "raw_intensity_imax_nA": raw_imax,
        "raw_intensity_hill_n": raw_hill,
        "filtered_intensity_k_mw_mm2": filt_k,
        "filtered_intensity_imax_nA": filt_imax,
        "filtered_intensity_hill_n": filt_hill,
        "soma_v_min_mV": float(np.min(trace_1s["soma_v"])),
        "soma_v_max_mV": float(np.max(trace_1s["soma_v"])),
        **area_summary(rows),
        **field,
    }
    return summary, intensity_rows, trace_1s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang_4state_no_truncate_rerun"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
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

    summaries = []
    all_intensity_rows = []
    traces = {}
    for mode in ["full", "apical_0p4mm2"]:
        summary, intensity_rows, trace = run_mode(args, protocol, base_rows, mode)
        summaries.append(summary)
        all_intensity_rows.extend(intensity_rows)
        traces[mode] = trace

    with (args.out_dir / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(summaries)
    with (args.out_dir / "intensity.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_intensity_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_intensity_rows)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2), constrained_layout=True)
    for mode, trace in traces.items():
        axes[0].plot(trace["t"], -trace["filtered_i"], lw=1.2, label=mode)
    axes[0].set_xlim(-10, 1100)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("current (nA)")
    axes[0].set_title("4-state 22C, 9.2 mW/mm2, no truncate")
    axes[0].legend(fontsize=8)

    for summary in summaries:
        mode = summary["mode"]
        sub = [row for row in all_intensity_rows if row["mode"] == mode]
        axes[1].plot(
            [row["irradiance_mw_mm2"] for row in sub],
            [row["filtered_peak_nA"] for row in sub],
            marker="o",
            label=f"{mode}, K={summary['filtered_intensity_k_mw_mm2']:.2g}",
        )
    axes[1].axvline(0.84, color="#d62728", ls="--", lw=1, label="Wang K~0.84")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("peak current (nA)")
    axes[1].set_title("Intensity response")
    axes[1].legend(fontsize=7)

    x = np.arange(len(summaries))
    axes[2].bar(x - 0.18, [row["retained_area_um2"] for row in summaries], width=0.36, label="retained")
    axes[2].bar(x + 0.18, [row["illuminated_area_um2"] for row in summaries], width=0.36, label="illuminated")
    axes[2].set_xticks(x, [row["mode"] for row in summaries], rotation=15)
    axes[2].set_ylabel("membrane area (um2)")
    axes[2].set_title("Area sanity check")
    axes[2].legend(fontsize=8)
    fig.savefig(args.out_dir / "summary.png", dpi=220)
    plt.close(fig)

    for summary in summaries:
        print(summary)


if __name__ == "__main__":
    main()
