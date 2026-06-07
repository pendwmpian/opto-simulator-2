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
from neuron import h

from calibrate_wang2007_gbar import Wang2007Protocol, current_metrics
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_pyrho6_fig2abc import (
    calibrate_pyrho6_gbar,
    q10_rate_scale,
    run_trace,
)
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    add_localization_weights,
    fit_hill_k,
    load_mechanisms,
    prepare_template_rows,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang_pyrho6_fig2b_traces"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--q10", type=float, default=1.0)
    parser.add_argument("--gamma-scale", type=float, default=5.0)
    parser.add_argument("--phi-m-scale", type=float, default=0.75)
    parser.add_argument("--rs-mohm", type=float, default=10.0)
    parser.add_argument("--electrode-rs-mohm", type=float, default=10.0)
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
    parser.add_argument("--binary-iterations", type=int, default=10)
    parser.add_argument("--pulse-duration-ms", type=float, default=100.0)
    parser.add_argument("--tstop-ms", type=float, default=180.0)
    parser.add_argument("--saturating-irradiance-mw-mm2", type=float, default=100.0)
    parser.add_argument("--saturating-duration-ms", type=float, default=100.0)
    parser.add_argument("--use-original-biophysics", action="store_true", default=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    protocol = Wang2007Protocol()
    h.load_file("stdrun.hoc")
    template_cell = CellFromNetPyNE(args.cell_params)
    rows = prepare_template_rows(template_cell, protocol, args)
    rows = add_localization_weights(rows, args.localization, args.proximal_cutoff_um)

    rate_scale = q10_rate_scale(args.q10, args.reference_c, args.target_c)
    gbar, calibration_peak = calibrate_pyrho6_gbar(
        args.cell_params,
        rows,
        args,
        protocol,
        rate_scale,
        args.gamma_scale,
        args.phi_m_scale,
    )

    metric_rows = []
    traces = {}
    for irr in WANG_INTENSITY_VALUES_MW_MM2:
        trace = run_trace(
            args.cell_params,
            rows,
            gbar,
            irr,
            args.pulse_duration_ms,
            args.tstop_ms,
            args,
            protocol,
            rate_scale,
            args.gamma_scale,
            args.phi_m_scale,
        )
        raw_peak, raw_steady, raw_ttp, raw_tau = current_metrics(trace["t"], trace["clamp_i"], protocol)
        filt_peak, filt_steady, filt_ttp, filt_tau = current_metrics(trace["t"], trace["filtered_i"], protocol)
        metric_rows.append(
            {
                "label": args.label,
                "q10": args.q10,
                "gamma_scale": args.gamma_scale,
                "phi_m_scale": args.phi_m_scale,
                "effective_gamma": 0.00369 * args.gamma_scale,
                "effective_phi_m": 5.02e17 * args.phi_m_scale,
                "gbar_mS_cm2": gbar,
                "calibration_peak_nA": calibration_peak,
                "irradiance_mw_mm2": irr,
                "raw_peak_nA": raw_peak,
                "raw_steady_nA": raw_steady,
                "raw_ttp_ms": raw_ttp,
                "raw_tau_ms": raw_tau,
                "filtered_peak_nA": filt_peak,
                "filtered_steady_nA": filt_steady,
                "filtered_ttp_ms": filt_ttp,
                "filtered_tau_ms": filt_tau,
            }
        )
        traces[irr] = trace

    intensity_k, intensity_imax, intensity_hill = fit_hill_k(
        [row["irradiance_mw_mm2"] for row in metric_rows],
        [row["filtered_peak_nA"] for row in metric_rows],
    )
    with (args.out_dir / "fig2b_100ms_intensity_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(metric_rows[0].keys()))
        writer.writeheader()
        writer.writerows(metric_rows)
    with (args.out_dir / "fig2b_100ms_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "label",
                "q10",
                "gamma_scale",
                "phi_m_scale",
                "gbar_mS_cm2",
                "intensity_k_mw_mm2",
                "intensity_imax_nA",
                "intensity_hill_n",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "label": args.label,
                "q10": args.q10,
                "gamma_scale": args.gamma_scale,
                "phi_m_scale": args.phi_m_scale,
                "gbar_mS_cm2": gbar,
                "intensity_k_mw_mm2": intensity_k,
                "intensity_imax_nA": intensity_imax,
                "intensity_hill_n": intensity_hill,
            }
        )

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for idx, irr in enumerate(WANG_INTENSITY_VALUES_MW_MM2):
        color = cmap(idx / (len(WANG_INTENSITY_VALUES_MW_MM2) - 1))
        trace = traces[irr]
        axes[0].plot(trace["t"], -trace["filtered_i"], color=color, lw=1.3, label=f"{irr:g}")
    axes[0].axhline(0.0, color="#333333", lw=0.6)
    axes[0].plot([0.0, args.pulse_duration_ms], [0.045, 0.045], color="#111111", lw=2.0)
    axes[0].set_xlim(-5.0, args.tstop_ms)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("current (nA)")
    axes[0].set_title("PyRhO6 Fig. 2B-like 100 ms traces")
    axes[0].legend(title="mW/mm2", fontsize=7, title_fontsize=8)

    axes[1].plot(
        [row["irradiance_mw_mm2"] for row in metric_rows],
        [row["filtered_peak_nA"] for row in metric_rows],
        marker="o",
        label=f"peak K={intensity_k:.2g}",
    )
    axes[1].plot(
        [row["irradiance_mw_mm2"] for row in metric_rows],
        [row["filtered_steady_nA"] for row in metric_rows],
        marker="s",
        label="end/steady readout",
    )
    axes[1].axvline(0.84, color="#d62728", ls="--", lw=1, label="Wang K~0.84")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("current (nA)")
    axes[1].set_title("Intensity response")
    axes[1].legend(fontsize=8)
    fig.savefig(args.out_dir / "fig2b_100ms_traces.png", dpi=220)
    plt.close(fig)

    for row in metric_rows:
        print(row)


if __name__ == "__main__":
    main()
