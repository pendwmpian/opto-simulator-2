from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from calibrate_wang2007_gbar import Wang2007Protocol, current_metrics, half_max_x
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    add_localization_weights,
    calibrate_gbar,
    fit_hill_k,
    load_mechanisms,
    prepare_template_rows,
    simulate_mod_current,
    williams_q10_scales,
)
from run_chr2_single_cell_sweep import CellFromNetPyNE
from neuron import h


def first_order_lowpass(t: np.ndarray, x: np.ndarray, tau_ms: float) -> np.ndarray:
    if tau_ms <= 0:
        return x.copy()
    y = np.empty_like(x)
    y[0] = x[0]
    for i in range(1, len(x)):
        dt = t[i] - t[i - 1]
        alpha = dt / (tau_ms + dt)
        y[i] = y[i - 1] + alpha * (x[i] - y[i - 1])
    return y


def apply_electrode_amplifier_readout(t: np.ndarray, current: np.ndarray, args) -> np.ndarray:
    tau_ms = args.electrode_rs_mohm * args.pipette_capacitance_pf * 1.0e-3
    y = first_order_lowpass(t, current, tau_ms)
    for _ in range(max(0, args.amplifier_filter_order - 1)):
        y = first_order_lowpass(t, y, args.amplifier_filter_tau_ms)
    return y


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
    filtered_i = apply_electrode_amplifier_readout(t, clamp_i, args)
    return {"t": t, "clamp_i": clamp_i, "filtered_i": filtered_i, "soma_v": soma_v}


def early_metrics(t: np.ndarray, current: np.ndarray, peak: float) -> tuple[float, float]:
    mask = (t >= 0.0) & (t <= 2.0)
    early = float(np.max(current[mask])) if np.any(mask) else float("nan")
    return early, early / peak if peak > 0.0 else float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_electrode_amplifier_model"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--rs-mohm", type=float, default=5.0)
    parser.add_argument("--electrode-rs-mohm", type=float, default=5.0)
    parser.add_argument("--pipette-capacitance-pf", type=float, default=100.0)
    parser.add_argument("--amplifier-filter-tau-ms", type=float, default=0.25)
    parser.add_argument("--amplifier-filter-order", type=int, default=4)
    parser.add_argument("--irradiance-scale", type=float, default=1.0)
    parser.add_argument("--soma-entry-depth-um", type=float, default=100.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--illumination-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=1.3)
    parser.add_argument("--localization", default="uniform")
    parser.add_argument("--proximal-cutoff-um", type=float, default=150.0)
    parser.add_argument("--reference-c", type=float, default=37.0)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--dt-ms", type=float, default=0.05)
    parser.add_argument("--tstop-ms", type=float, default=180.0)
    parser.add_argument("--pulse-duration-ms", type=float, default=100.0)
    parser.add_argument("--binary-iterations", type=int, default=12)
    parser.add_argument("--saturating-irradiance-mw-mm2", type=float, default=100.0)
    parser.add_argument("--saturating-duration-ms", type=float, default=100.0)
    parser.add_argument("--c2-initial", type=float, default=0.0)
    parser.add_argument("--e12-scale", type=float, default=1.0)
    parser.add_argument("--e21-scale", type=float, default=1.0)
    parser.add_argument("--gd2-scale", type=float, default=1.0)
    parser.add_argument("--gamma-scale", type=float, default=1.0)
    parser.add_argument("--use-original-biophysics", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    protocol = Wang2007Protocol()
    h.load_file("stdrun.hoc")
    template_cell = CellFromNetPyNE(args.cell_params)
    rows = prepare_template_rows(template_cell, protocol, args)
    rows = add_localization_weights(rows, args.localization, args.proximal_cutoff_um)
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

    trace_9p2 = run_trace(args.cell_params, rows, gbar, 9.2, 1000.0, 1200.0, args, protocol, q10_scales)
    raw_peak, _, raw_ttp, raw_tau = current_metrics(trace_9p2["t"], trace_9p2["clamp_i"], protocol)
    filt_peak, _, filt_ttp, filt_tau = current_metrics(trace_9p2["t"], trace_9p2["filtered_i"], protocol)
    raw_early, raw_early_frac = early_metrics(trace_9p2["t"], trace_9p2["clamp_i"], raw_peak)
    filt_early, filt_early_frac = early_metrics(trace_9p2["t"], trace_9p2["filtered_i"], filt_peak)

    intensity_peaks = []
    intensity_ttps = []
    raw_intensity_peaks = []
    for irr in WANG_INTENSITY_VALUES_MW_MM2:
        trace = run_trace(args.cell_params, rows, gbar, irr, 100.0, 180.0, args, protocol, q10_scales)
        raw_peak_i, _, _, _ = current_metrics(trace["t"], trace["clamp_i"], protocol)
        peak_i, _, ttp_i, _ = current_metrics(trace["t"], trace["filtered_i"], protocol)
        raw_intensity_peaks.append(float(raw_peak_i))
        intensity_peaks.append(float(peak_i))
        intensity_ttps.append(float(ttp_i))
    intensity_k, intensity_imax, intensity_hill_n = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, intensity_peaks)
    raw_intensity_k, raw_intensity_imax, raw_intensity_hill_n = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, raw_intensity_peaks)

    duration_values = [1, 2, 3, 4, 5, 8, 10, 20, 50, 100]
    duration_peaks = []
    for dur in duration_values:
        trace = run_trace(args.cell_params, rows, gbar, 9.2, float(dur), max(120.0, float(dur) + 60.0), args, protocol, q10_scales)
        peak_d, _, _, _ = current_metrics(trace["t"], trace["filtered_i"], protocol)
        duration_peaks.append(float(peak_d))
    duration_k = half_max_x([float(v) for v in duration_values], duration_peaks)

    row = {
        "label": args.label,
        "gbar_mS_cm2": gbar,
        "rs_mohm": args.rs_mohm,
        "electrode_rs_mohm": args.electrode_rs_mohm,
        "pipette_capacitance_pf": args.pipette_capacitance_pf,
        "electrode_tau_ms": args.electrode_rs_mohm * args.pipette_capacitance_pf * 1.0e-3,
        "amplifier_filter_tau_ms": args.amplifier_filter_tau_ms,
        "amplifier_filter_order": args.amplifier_filter_order,
        "raw_peak_9p2_nA": raw_peak,
        "raw_ttp_9p2_ms": raw_ttp,
        "raw_tau_9p2_ms": raw_tau,
        "raw_early_peak_fraction": raw_early_frac,
        "filtered_peak_9p2_nA": filt_peak,
        "filtered_ttp_9p2_ms": filt_ttp,
        "filtered_tau_9p2_ms": filt_tau,
        "filtered_early_peak_fraction": filt_early_frac,
        "raw_intensity_k_mw_mm2": raw_intensity_k,
        "raw_intensity_imax_nA": raw_intensity_imax,
        "raw_intensity_hill_n": raw_intensity_hill_n,
        "filtered_intensity_k_mw_mm2": intensity_k,
        "filtered_intensity_imax_nA": intensity_imax,
        "filtered_intensity_hill_n": intensity_hill_n,
        "filtered_duration_k_ms": duration_k,
        "soma_v_min_mV": float(np.min(trace_9p2["soma_v"])),
        "soma_v_max_mV": float(np.max(trace_9p2["soma_v"])),
        "soma_v_mean_mV": float(np.mean(trace_9p2["soma_v"])),
        "filtered_ttp_0p07_ms": intensity_ttps[0],
        "filtered_ttp_0p58_ms": intensity_ttps[3],
        "filtered_ttp_2p3_ms": intensity_ttps[5],
        "filtered_ttp_9p2_ms": intensity_ttps[6],
    }
    with (args.out_dir / "electrode_amplifier_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    axes[0].plot(trace_9p2["t"], trace_9p2["clamp_i"], color="#777777", lw=1.1, label="raw SEClamp")
    axes[0].plot(trace_9p2["t"], trace_9p2["filtered_i"], color="#1f77b4", lw=1.4, label="electrode/amplifier readout")
    axes[0].set_xlim(-1, 60)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("current (nA)")
    axes[0].set_title("9.2 mW/mm2, 1 s")
    axes[0].legend(fontsize=8)
    axes[1].plot(WANG_INTENSITY_VALUES_MW_MM2, raw_intensity_peaks, marker="o", color="#777777", label="raw")
    axes[1].plot(WANG_INTENSITY_VALUES_MW_MM2, intensity_peaks, marker="o", color="#1f77b4", label="filtered")
    axes[1].axvline(0.84, color="#d62728", ls="--", lw=1)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("peak current (nA)")
    axes[1].set_title(f"K={intensity_k:.2f} mW/mm2")
    axes[1].legend(fontsize=8)
    axes[2].plot(duration_values, duration_peaks, marker="o", color="#1f77b4")
    axes[2].axvline(3.2, color="#d62728", ls="--", lw=1)
    axes[2].set_xscale("log")
    axes[2].set_xlabel("duration (ms)")
    axes[2].set_ylabel("peak current (nA)")
    axes[2].set_title(f"duration K={duration_k:.2f} ms")
    fig.savefig(args.out_dir / "electrode_amplifier_summary.png", dpi=220)
    plt.close(fig)

    print(row)


if __name__ == "__main__":
    main()
