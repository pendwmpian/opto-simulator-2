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

from calibrate_wang2007_gbar import Wang2007Protocol, half_max_x
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_electrode_amplifier_model import apply_electrode_amplifier_readout, run_trace
from run_wang_seclamp_mod_suite import (
    add_localization_weights,
    load_mechanisms,
    prepare_template_rows,
    williams_q10_scales,
)


def parse_floats(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def short_pulse_metrics(t: np.ndarray, current: np.ndarray, pulse_duration_ms: float) -> dict[str, float]:
    peak_idx = int(np.argmax(current))
    peak = float(current[peak_idx])
    ttp = float(t[peak_idx])
    end_mask = (t >= pulse_duration_ms - 0.25) & (t <= pulse_duration_ms + 0.25)
    post_mask = (t >= pulse_duration_ms) & (t <= min(float(t[-1]), pulse_duration_ms + 50.0))
    charge_mask = (t >= 0.0) & (t <= min(float(t[-1]), pulse_duration_ms + 50.0))
    end_current = float(np.mean(current[end_mask])) if np.any(end_mask) else float("nan")
    post_peak = float(np.max(current[post_mask])) if np.any(post_mask) else float("nan")
    charge_nC = float(np.trapezoid(current[charge_mask], t[charge_mask]) / 1000.0) if np.any(charge_mask) else float("nan")
    return {
        "peak_current_nA": peak,
        "time_to_peak_ms": ttp,
        "current_at_light_off_nA": end_current,
        "post_light_peak_nA": post_peak,
        "charge_0_to_off_plus_50ms_nC": charge_nC,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang_fig2de_short_pulse"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--gbar-mS-cm2", type=float, default=0.10090350448414476)
    parser.add_argument("--irradiance-mw-mm2", type=float, default=9.2)
    parser.add_argument("--durations-ms", default="1,2,3,4,5,6,7,8,10,20,30,40,50,60,70,80,90,100")
    parser.add_argument("--trace-durations-ms", default="1,2,3,4,5,6,7,8,10,20,30,40,50,60,70,80,90,100")
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
    parser.add_argument("--dt-ms", type=float, default=0.025)
    parser.add_argument("--tstop-ms", type=float, default=120.0)
    parser.add_argument("--use-original-biophysics", action="store_true", default=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    protocol = Wang2007Protocol()
    h.load_file("stdrun.hoc")
    template_cell = CellFromNetPyNE(args.cell_params)
    rows = prepare_template_rows(template_cell, protocol, args)
    rows = add_localization_weights(rows, args.localization, args.proximal_cutoff_um)
    q10_scales = williams_q10_scales(args.reference_c, args.target_c)

    durations = parse_floats(args.durations_ms)
    trace_durations = set(parse_floats(args.trace_durations_ms))
    metric_rows = []
    traces = {}
    for duration in durations:
        tstop = max(args.tstop_ms, duration + 80.0)
        trace = run_trace(
            args.cell_params,
            rows,
            args.gbar_mS_cm2,
            args.irradiance_mw_mm2,
            duration,
            tstop,
            args,
            protocol,
            q10_scales,
        )
        raw_metrics = short_pulse_metrics(trace["t"], trace["clamp_i"], duration)
        filtered_metrics = short_pulse_metrics(trace["t"], trace["filtered_i"], duration)
        metric_rows.append(
            {
                "label": args.label,
                "duration_ms": duration,
                "irradiance_mw_mm2": args.irradiance_mw_mm2,
                "gbar_mS_cm2": args.gbar_mS_cm2,
                "raw_peak_current_nA": raw_metrics["peak_current_nA"],
                "raw_time_to_peak_ms": raw_metrics["time_to_peak_ms"],
                "raw_charge_nC": raw_metrics["charge_0_to_off_plus_50ms_nC"],
                "filtered_peak_current_nA": filtered_metrics["peak_current_nA"],
                "filtered_time_to_peak_ms": filtered_metrics["time_to_peak_ms"],
                "filtered_current_at_light_off_nA": filtered_metrics["current_at_light_off_nA"],
                "filtered_post_light_peak_nA": filtered_metrics["post_light_peak_nA"],
                "filtered_charge_nC": filtered_metrics["charge_0_to_off_plus_50ms_nC"],
                "soma_v_min_mV": float(np.min(trace["soma_v"])),
                "soma_v_max_mV": float(np.max(trace["soma_v"])),
            }
        )
        if duration in trace_durations:
            traces[duration] = trace

    peak_max = max(row["filtered_peak_current_nA"] for row in metric_rows)
    charge_max = max(row["filtered_charge_nC"] for row in metric_rows)
    for row in metric_rows:
        row["filtered_peak_norm"] = row["filtered_peak_current_nA"] / peak_max if peak_max else float("nan")
        row["filtered_charge_norm"] = row["filtered_charge_nC"] / charge_max if charge_max else float("nan")

    duration_k_peak = half_max_x([row["duration_ms"] for row in metric_rows], [row["filtered_peak_current_nA"] for row in metric_rows])
    duration_k_charge = half_max_x([row["duration_ms"] for row in metric_rows], [row["filtered_charge_nC"] for row in metric_rows])

    with (args.out_dir / "fig2de_short_pulse_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(metric_rows[0].keys()))
        writer.writeheader()
        writer.writerows(metric_rows)

    summary = {
        "label": args.label,
        "gbar_mS_cm2": args.gbar_mS_cm2,
        "irradiance_mw_mm2": args.irradiance_mw_mm2,
        "duration_k_peak_ms": duration_k_peak,
        "duration_k_charge_ms": duration_k_charge,
        "rs_mohm": args.rs_mohm,
        "electrode_rs_mohm": args.electrode_rs_mohm,
        "pipette_capacitance_pf": args.pipette_capacitance_pf,
        "amplifier_filter_tau_ms": args.amplifier_filter_tau_ms,
        "amplifier_filter_order": args.amplifier_filter_order,
        "vitro_mu_eff_mm_inv": args.vitro_mu_eff_mm_inv,
    }
    with (args.out_dir / "fig2de_short_pulse_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary.keys()))
        writer.writeheader()
        writer.writerow(summary)

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    d_traces = [(duration, trace) for duration, trace in sorted(traces.items()) if duration <= 8.0]
    for idx, (duration, trace) in enumerate(d_traces):
        color = cmap(idx / max(1, len(d_traces) - 1))
        axes[0].plot(trace["t"], -trace["filtered_i"], lw=1.25, color=color, label=f"{duration:g} ms")
        axes[0].plot([0.0, duration], [0.035 + idx * 0.012, 0.035 + idx * 0.012], color=color, lw=2.0)
    axes[0].axhline(0.0, color="#333333", lw=0.6)
    axes[0].set_xlim(-1.0, 40.0)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("current (nA)")
    axes[0].set_title("Fig. 2D-like 1-8 ms flashes")
    axes[0].legend(fontsize=7, ncol=2, loc="lower right")

    e_traces = [(duration, trace) for duration, trace in sorted(traces.items()) if duration >= 10.0]
    e_colors = plt.get_cmap("plasma")
    for idx, (duration, trace) in enumerate(e_traces):
        color = e_colors(idx / max(1, len(e_traces) - 1))
        axes[1].plot(trace["t"], -trace["filtered_i"], lw=1.1, color=color, label=f"{duration:g} ms")
    axes[1].axhline(0.0, color="#333333", lw=0.6)
    axes[1].set_xlim(-2.0, 150.0)
    axes[1].set_xlabel("time (ms)")
    axes[1].set_ylabel("current (nA)")
    axes[1].set_title("Fig. 2E-like 10-100 ms flashes")
    axes[1].legend(fontsize=6, ncol=2, loc="lower right")

    xs = [row["duration_ms"] for row in metric_rows]
    axes[2].plot(xs, [row["filtered_peak_norm"] for row in metric_rows], marker="o", color="#1f77b4", label="peak")
    axes[2].plot(xs, [row["filtered_charge_norm"] for row in metric_rows], marker="s", color="#2ca02c", label="charge")
    if duration_k_peak is not None:
        axes[2].axvline(duration_k_peak, color="#1f77b4", ls="--", lw=1)
    if duration_k_charge is not None:
        axes[2].axvline(duration_k_charge, color="#2ca02c", ls=":", lw=1)
    axes[2].set_xscale("log")
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].set_xlabel("light duration (ms)")
    axes[2].set_ylabel("normalized response")
    axes[2].set_title("Duration-response summary")
    axes[2].legend(fontsize=8)
    fig.savefig(args.out_dir / "fig2de_short_pulse.png", dpi=220)
    plt.close(fig)

    for row in metric_rows:
        print(row)
    print(summary)


if __name__ == "__main__":
    main()
