from __future__ import annotations

import argparse
import csv
import itertools
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from neuron import h

from calibrate_wang2007_gbar import Wang2007Protocol, current_metrics, half_max_x
from run_wang_electrode_amplifier_model import apply_electrode_amplifier_readout, run_trace
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


def parse_values(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def score(row: dict[str, float]) -> float:
    targets = {
        "filtered_peak_9p2_nA": (0.557, 0.153),
        "filtered_ttp_9p2_ms": (12.0, 3.0),
        "filtered_tau_9p2_ms": (48.0, 12.0),
        "filtered_intensity_k_mw_mm2": (0.84, 0.2),
    }
    if np.isfinite(row.get("filtered_duration_k_ms", float("nan"))):
        targets["filtered_duration_k_ms"] = (3.2, 0.8)
    total = 0.0
    for key, (target, scale) in targets.items():
        value = row[key]
        if value is None or not np.isfinite(value):
            total += 100.0
        else:
            total += ((value - target) / scale) ** 2
    return float(total)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_desensitization_sweep"))
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
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--localization", default="uniform")
    parser.add_argument("--proximal-cutoff-um", type=float, default=150.0)
    parser.add_argument("--reference-c", type=float, default=37.0)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--tstop-ms", type=float, default=180.0)
    parser.add_argument("--pulse-duration-ms", type=float, default=100.0)
    parser.add_argument("--binary-iterations", type=int, default=9)
    parser.add_argument("--saturating-irradiance-mw-mm2", type=float, default=100.0)
    parser.add_argument("--saturating-duration-ms", type=float, default=100.0)
    parser.add_argument("--c2-initial", type=float, default=0.0)
    parser.add_argument("--e12-scales", default="0.5,1.0")
    parser.add_argument("--e21-scales", default="1.0,2.0")
    parser.add_argument("--gd2-scales", default="0.5,1.0")
    parser.add_argument("--gamma-scales", default="1.0,2.0")
    parser.add_argument("--skip-duration", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    protocol = Wang2007Protocol()
    h.load_file("stdrun.hoc")
    template_cell = CellFromNetPyNE(args.cell_params)
    rows = prepare_template_rows(template_cell, protocol, args)
    rows = add_localization_weights(rows, args.localization, args.proximal_cutoff_um)
    q10_scales = williams_q10_scales(args.reference_c, args.target_c)

    results = []
    grid = itertools.product(
        parse_values(args.e12_scales),
        parse_values(args.e21_scales),
        parse_values(args.gd2_scales),
        parse_values(args.gamma_scales),
    )
    for e12_scale, e21_scale, gd2_scale, gamma_scale in grid:
        args.e12_scale = e12_scale
        args.e21_scale = e21_scale
        args.gd2_scale = gd2_scale
        args.gamma_scale = gamma_scale
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

        trace_9p2 = run_trace(args.cell_params, rows, gbar, 9.2, 1000.0, 1200.0, args, protocol)
        filt_peak, _, filt_ttp, filt_tau = current_metrics(trace_9p2["t"], trace_9p2["filtered_i"], protocol)

        intensity_peaks = []
        for irr in WANG_INTENSITY_VALUES_MW_MM2:
            trace = run_trace(args.cell_params, rows, gbar, irr, 100.0, 180.0, args, protocol)
            peak_i, _, _, _ = current_metrics(trace["t"], trace["filtered_i"], protocol)
            intensity_peaks.append(float(peak_i))
        intensity_k, intensity_imax, intensity_hill_n = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, intensity_peaks)

        duration_k = float("nan")
        if not args.skip_duration:
            duration_values = [1, 2, 3, 4, 5, 8, 10, 20, 50, 100]
            duration_peaks = []
            for dur in duration_values:
                trace = run_trace(args.cell_params, rows, gbar, 9.2, float(dur), max(120.0, float(dur) + 60.0), args, protocol)
                peak_d, _, _, _ = current_metrics(trace["t"], trace["filtered_i"], protocol)
                duration_peaks.append(float(peak_d))
            duration_k = half_max_x([float(v) for v in duration_values], duration_peaks)

        row = {
            "label": args.label,
            "e12_scale": e12_scale,
            "e21_scale": e21_scale,
            "gd2_scale": gd2_scale,
            "gamma_scale": gamma_scale,
            "gbar_mS_cm2": float(gbar),
            "filtered_peak_9p2_nA": float(filt_peak),
            "filtered_ttp_9p2_ms": float(filt_ttp),
            "filtered_tau_9p2_ms": float(filt_tau) if filt_tau is not None else float("nan"),
            "filtered_intensity_k_mw_mm2": float(intensity_k) if intensity_k is not None else float("nan"),
            "filtered_intensity_imax_nA": float(intensity_imax),
            "filtered_intensity_hill_n": float(intensity_hill_n),
            "filtered_duration_k_ms": float(duration_k) if duration_k is not None else float("nan"),
        }
        row["score"] = score(row)
        results.append(row)
        print(row)

    results.sort(key=lambda item: item["score"])
    with (args.out_dir / "desensitization_sweep.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    top = results[: min(8, len(results))]
    labels = [f"e12 {r['e12_scale']}, e21 {r['e21_scale']}, gd2 {r['gd2_scale']}, g {r['gamma_scale']}" for r in top]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), constrained_layout=True)
    metrics = [
        ("filtered_ttp_9p2_ms", 12.0, "time-to-peak (ms)"),
        ("filtered_tau_9p2_ms", 48.0, "tau (ms)"),
        ("filtered_intensity_k_mw_mm2", 0.84, "intensity K"),
        ("filtered_duration_k_ms", 3.2, "duration K"),
    ]
    for ax, (key, target, title) in zip(axes, metrics):
        ax.bar(range(len(top)), [r[key] for r in top], color="#4c78a8")
        ax.axhline(target, color="#d62728", ls="--", lw=1)
        ax.set_title(title)
        ax.set_xticks(range(len(top)))
        ax.set_xticklabels(labels, rotation=80, ha="right", fontsize=7)
    fig.savefig(args.out_dir / "desensitization_sweep_top.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
