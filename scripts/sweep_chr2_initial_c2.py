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

from calibrate_wang2007_gbar import Wang2007Protocol, current_metrics
from diagnose_wang_chr2_intensity_kinetics import simulate_diagnostic
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    add_localization_weights,
    calibrate_gbar,
    fit_hill_k,
    load_mechanisms,
    prepare_template_rows,
    williams_q10_scales,
)
from run_chr2_single_cell_sweep import CellFromNetPyNE
from neuron import h


def early_transient_metrics(t: np.ndarray, current: np.ndarray, peak: float) -> dict[str, float]:
    early_mask = (t >= 0.0) & (t <= 2.0)
    mid_mask = (t >= 2.0) & (t <= 12.0)
    if not np.any(early_mask) or not np.any(mid_mask) or peak <= 0.0:
        return {"early_peak_nA": float("nan"), "early_peak_fraction": float("nan"), "early_step_score": float("nan")}
    early_peak = float(np.max(current[early_mask]))
    mid_peak = float(np.max(current[mid_mask]))
    return {
        "early_peak_nA": early_peak,
        "early_peak_fraction": early_peak / peak,
        "early_step_score": early_peak / mid_peak if mid_peak > 0.0 else float("nan"),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_c2_initial_sweep"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--rs-mohm", type=float, default=5.0)
    parser.add_argument("--irradiance-scale", type=float, default=1.0)
    parser.add_argument("--soma-entry-depth-um", type=float, default=100.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--illumination-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--localization", default="uniform")
    parser.add_argument("--proximal-cutoff-um", type=float, default=150.0)
    parser.add_argument("--reference-c", type=float, default=37.0)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--tstop-ms", type=float, default=180.0)
    parser.add_argument("--pulse-duration-ms", type=float, default=100.0)
    parser.add_argument("--binary-iterations", type=int, default=12)
    parser.add_argument("--saturating-irradiance-mw-mm2", type=float, default=100.0)
    parser.add_argument("--saturating-duration-ms", type=float, default=100.0)
    parser.add_argument("--c2-values", default="0,0.05,0.1,0.2,0.3,0.5")
    parser.add_argument("--trace-irradiance-mw-mm2", type=float, default=9.2)
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
    )

    c2_values = [float(raw.strip()) for raw in args.c2_values.split(",") if raw.strip()]
    summary = []
    trace_by_c2 = {}
    for c2_init in c2_values:
        args.c2_initial = c2_init
        intensity_peaks = []
        intensity_ttps = []
        for irradiance in WANG_INTENSITY_VALUES_MW_MM2:
            trace = simulate_diagnostic(args.cell_params, rows, gbar, irradiance, args, protocol)
            peak, _, ttp, tau = current_metrics(trace["t"], trace["clamp_i"], protocol)
            intensity_peaks.append(float(peak))
            intensity_ttps.append(float(ttp))
            if abs(irradiance - args.trace_irradiance_mw_mm2) < 1.0e-9:
                trace_by_c2[c2_init] = trace
                trace_peak = float(peak)
                trace_ttp = float(ttp)
                trace_tau = tau
                trace_early = early_transient_metrics(trace["t"], trace["clamp_i"], trace_peak)
        intensity_k, intensity_imax, intensity_hill_n = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, intensity_peaks)
        summary.append(
            {
                "label": args.label,
                "c1_initial": 1.0 - c2_init,
                "c2_initial": c2_init,
                "gbar_mS_cm2": gbar,
                "trace_irradiance_mw_mm2": args.trace_irradiance_mw_mm2,
                "trace_peak_nA": trace_peak,
                "trace_time_to_peak_ms": trace_ttp,
                "trace_tau_ms": trace_tau,
                **trace_early,
                "intensity_k_mw_mm2": intensity_k,
                "intensity_imax_nA": intensity_imax,
                "intensity_hill_n": intensity_hill_n,
                "ttp_0p07_ms": intensity_ttps[0],
                "ttp_0p58_ms": intensity_ttps[3],
                "ttp_2p3_ms": intensity_ttps[5],
                "ttp_9p2_ms": intensity_ttps[6],
            }
        )

    with (args.out_dir / "c2_initial_sweep.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    for c2_init, trace in trace_by_c2.items():
        axes[0].plot(trace["t"], trace["clamp_i"], lw=1.2, label=f"C2={c2_init:g}")
        axes[1].plot(trace["t"], trace["open"], lw=1.2)
    axes[0].set_xlim(-1, 60)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("SEClamp current (nA)")
    axes[0].set_title(f"{args.trace_irradiance_mw_mm2:g} mW/mm2 current")
    axes[0].legend(fontsize=8)
    axes[1].set_xlim(-1, 60)
    axes[1].set_xlabel("time (ms)")
    axes[1].set_ylabel("weighted open")
    axes[1].set_title("ChR2 open state")
    axes[2].plot([row["c2_initial"] for row in summary], [row["early_peak_fraction"] for row in summary], marker="o", label="early peak / total peak")
    axes[2].plot([row["c2_initial"] for row in summary], [row["trace_time_to_peak_ms"] for row in summary], marker="o", label="time to peak (ms)")
    axes[2].set_xlabel("initial C2 fraction")
    axes[2].set_title("early transient suppression")
    axes[2].legend(fontsize=8)
    fig.savefig(args.out_dir / "c2_initial_sweep.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
    ax.plot([row["c2_initial"] for row in summary], [row["intensity_k_mw_mm2"] for row in summary], marker="o", label="intensity K")
    ax.axhline(0.84, color="#d62728", ls="--", lw=1, label="Wang K")
    ax.set_xlabel("initial C2 fraction")
    ax.set_ylabel("mW/mm2")
    ax.legend()
    fig.savefig(args.out_dir / "c2_initial_intensity_k.png", dpi=220)
    plt.close(fig)

    for row in summary:
        print(row)


if __name__ == "__main__":
    main()
