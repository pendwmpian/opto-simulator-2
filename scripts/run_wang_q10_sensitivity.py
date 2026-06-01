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

from calibrate_wang2007_gbar import Wang2007Protocol, current_metrics, half_max_x
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_electrode_amplifier_model import apply_electrode_amplifier_readout, run_trace
from run_wang_fig2de_short_pulse import short_pulse_metrics
from run_wang_seclamp_mod_suite import (
    add_localization_weights,
    calibrate_gbar,
    generic_q10_scales,
    load_mechanisms,
    prepare_template_rows,
    williams_q10_scales,
)


def q10_condition_scales(name: str, reference_c: float, target_c: float) -> dict[str, float]:
    if name == "none":
        return {key: 1.0 for key in ["k1", "k2", "gd1", "gd2", "e12", "e21", "gr"]}
    if name == "williams":
        return williams_q10_scales(reference_c, target_c)
    if name.startswith("generic"):
        q10 = float(name.replace("generic", ""))
        return generic_q10_scales(q10, reference_c, target_c)
    raise ValueError(f"Unknown q10 condition: {name}")


def parse_items(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def parse_floats(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang_q10_sensitivity"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--q10-conditions", default="none,generic1.5,generic2.0,generic2.5,williams")
    parser.add_argument("--durations-ms", default="1,2,3,4,5,6,7,8,10,20,30,40,50,60,70,80,90,100")
    parser.add_argument("--irradiance-mw-mm2", type=float, default=9.2)
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
    parser.add_argument("--binary-iterations", type=int, default=8)
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

    summaries = []
    duration_rows = []
    duration_values = parse_floats(args.durations_ms)
    for condition in parse_items(args.q10_conditions):
        q10_scales = q10_condition_scales(condition, args.reference_c, args.target_c)
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

        trace_1s = run_trace(
            args.cell_params,
            rows,
            gbar,
            args.irradiance_mw_mm2,
            1000.0,
            1200.0,
            args,
            protocol,
            q10_scales,
        )
        raw_peak, raw_steady, raw_ttp, raw_tau = current_metrics(trace_1s["t"], trace_1s["clamp_i"], protocol)
        filtered_peak, filtered_steady, filtered_ttp, filtered_tau = current_metrics(trace_1s["t"], trace_1s["filtered_i"], protocol)

        peaks = []
        charges = []
        for duration in duration_values:
            trace = run_trace(
                args.cell_params,
                rows,
                gbar,
                args.irradiance_mw_mm2,
                duration,
                max(120.0, duration + 80.0),
                args,
                protocol,
                q10_scales,
            )
            metrics = short_pulse_metrics(trace["t"], trace["filtered_i"], duration)
            peaks.append(metrics["peak_current_nA"])
            charges.append(metrics["charge_0_to_off_plus_50ms_nC"])
            duration_rows.append(
                {
                    "q10_condition": condition,
                    "duration_ms": duration,
                    "gbar_mS_cm2": gbar,
                    "filtered_peak_current_nA": metrics["peak_current_nA"],
                    "filtered_time_to_peak_ms": metrics["time_to_peak_ms"],
                    "filtered_current_at_light_off_nA": metrics["current_at_light_off_nA"],
                    "filtered_charge_nC": metrics["charge_0_to_off_plus_50ms_nC"],
                }
            )
        peak_max = max(peaks)
        charge_max = max(charges)
        for row in duration_rows:
            if row["q10_condition"] != condition:
                continue
            row["filtered_peak_norm"] = row["filtered_peak_current_nA"] / peak_max if peak_max else float("nan")
            row["filtered_charge_norm"] = row["filtered_charge_nC"] / charge_max if charge_max else float("nan")

        summaries.append(
            {
                "q10_condition": condition,
                "gbar_mS_cm2": gbar,
                "k1_scale": q10_scales["k1"],
                "k2_scale": q10_scales["k2"],
                "gd1_scale": q10_scales["gd1"],
                "gd2_scale": q10_scales["gd2"],
                "e12_scale": q10_scales["e12"],
                "e21_scale": q10_scales["e21"],
                "gr_scale": q10_scales["gr"],
                "raw_peak_1s_nA": raw_peak,
                "raw_steady_1s_nA": raw_steady,
                "raw_steady_peak_ratio": raw_steady / raw_peak if raw_peak else float("nan"),
                "raw_ttp_1s_ms": raw_ttp,
                "raw_tau_1s_ms": raw_tau,
                "filtered_peak_1s_nA": filtered_peak,
                "filtered_steady_1s_nA": filtered_steady,
                "filtered_steady_peak_ratio": filtered_steady / filtered_peak if filtered_peak else float("nan"),
                "filtered_ttp_1s_ms": filtered_ttp,
                "filtered_tau_1s_ms": filtered_tau,
                "duration_k_peak_ms": half_max_x(duration_values, peaks),
                "duration_k_charge_ms": half_max_x(duration_values, charges),
            }
        )

    with (args.out_dir / "q10_sensitivity_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(summaries)
    with (args.out_dir / "q10_sensitivity_duration.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(duration_rows[0].keys()))
        writer.writeheader()
        writer.writerows(duration_rows)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2), constrained_layout=True)
    for summary in summaries:
        condition = summary["q10_condition"]
        sub = [row for row in duration_rows if row["q10_condition"] == condition]
        axes[0].plot([row["duration_ms"] for row in sub], [row["filtered_peak_norm"] for row in sub], marker="o", label=condition)
        axes[1].plot([row["duration_ms"] for row in sub], [row["filtered_charge_norm"] for row in sub], marker="o", label=condition)
    for ax in axes[:2]:
        ax.set_xscale("log")
        ax.set_xlabel("duration (ms)")
        ax.set_ylim(-0.05, 1.05)
    axes[0].set_ylabel("normalized peak")
    axes[0].set_title("Peak duration response")
    axes[1].set_ylabel("normalized charge")
    axes[1].set_title("Charge duration response")
    axes[1].legend(fontsize=7)

    x = np.arange(len(summaries))
    axes[2].bar(x - 0.18, [row["filtered_peak_1s_nA"] for row in summaries], width=0.36, label="peak")
    axes[2].bar(x + 0.18, [row["filtered_steady_1s_nA"] for row in summaries], width=0.36, label="steady")
    axes[2].axhline(0.4, color="#d62728", ls="--", lw=1, label="Wang ~0.4 nA")
    axes[2].set_xticks(x, [row["q10_condition"] for row in summaries], rotation=25, ha="right")
    axes[2].set_ylabel("current (nA)")
    axes[2].set_title("1 s photocurrent")
    axes[2].legend(fontsize=8)
    fig.savefig(args.out_dir / "q10_sensitivity.png", dpi=220)
    plt.close(fig)

    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
