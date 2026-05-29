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
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    add_localization_weights,
    calibrate_gbar,
    install_chr2,
    load_mechanisms,
    prepare_template_rows,
    soma_section,
    williams_q10_scales,
)
from run_chr2_single_cell_sweep import CellFromNetPyNE


def weighted_state(rows, attr: str) -> float:
    numer = 0.0
    denom = 0.0
    for row in rows:
        weight = row["area_um2"] * max(float(row["seg"].chr2_4state.gbar), 0.0)
        numer += weight * float(getattr(row["seg"].chr2_4state, attr))
        denom += weight
    return numer / denom if denom > 0.0 else 0.0


def summed_chr2_current_nA(rows) -> float:
    total = 0.0
    for row in rows:
        # i is mA/cm2 and area is um2. 1 um2 = 1e-8 cm2, 1 mA = 1e6 nA.
        total += abs(float(row["seg"].chr2_4state.i)) * row["area_um2"] * 1.0e-2
    return total


def first_peak_time(t: np.ndarray, y: np.ndarray, max_t_ms: float = 120.0) -> float:
    mask = t <= max_t_ms
    if not np.any(mask):
        return float("nan")
    idx = int(np.argmax(y[mask]))
    return float(t[mask][idx])


def simulate_diagnostic(cell_params, rows_template, gbar, irradiance, args, protocol):
    h.load_file("stdrun.hoc")
    h.celsius = args.target_c
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=getattr(args, "use_original_biophysics", False))
    q10_scales = williams_q10_scales(args.reference_c, args.target_c)
    rows = install_chr2(cell, rows_template, gbar, args.irradiance_scale, q10_scales, True)
    c2_init = args.c2_initial
    c1_init = max(0.0, 1.0 - c2_init)
    for row in rows:
        mech = row["seg"].chr2_4state
        mech.c1_init = c1_init
        mech.c2_init = c2_init
        mech.o1_init = 0.0
        mech.o2_init = 0.0
        mech.fit_e12_scale = getattr(args, "e12_scale", 1.0)
        mech.fit_e21_scale = getattr(args, "e21_scale", 1.0)
        mech.fit_gd2_scale = getattr(args, "gd2_scale", 1.0)
        mech.fit_gamma_scale = getattr(args, "gamma_scale", 1.0)

    clamp = h.SEClamp(soma_section(cell)(0.5))
    clamp.dur1 = args.tstop_ms + 1.0
    clamp.amp1 = protocol.holding_voltage_mV
    clamp.rs = args.rs_mohm

    h.dt = args.dt_ms
    h.tstop = args.tstop_ms
    h.finitialize(protocol.holding_voltage_mV)

    times = []
    clamp_i = []
    local_i = []
    open_state = []
    o1_state = []
    o2_state = []
    c2_state = []
    p_state = []
    while h.t <= args.tstop_ms:
        light_on = 0.0 <= h.t < args.pulse_duration_ms
        scale = irradiance / protocol.irradiance_mw_mm2
        for row in rows:
            row["seg"].chr2_4state.irr = row["base_irradiance_mw_mm2"] * scale if light_on else 0.0
        times.append(float(h.t))
        clamp_i.append(abs(float(clamp.i)))
        local_i.append(summed_chr2_current_nA(rows))
        open_state.append(weighted_state(rows, "open"))
        o1_state.append(weighted_state(rows, "o1"))
        o2_state.append(weighted_state(rows, "o2"))
        c2_state.append(weighted_state(rows, "c2"))
        p_state.append(weighted_state(rows, "p"))
        h.fadvance()

    return {
        "t": np.asarray(times),
        "clamp_i": np.asarray(clamp_i),
        "local_i": np.asarray(local_i),
        "open": np.asarray(open_state),
        "o1": np.asarray(o1_state),
        "o2": np.asarray(o2_state),
        "c2": np.asarray(c2_state),
        "p": np.asarray(p_state),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_chr2_intensity_diagnostics"))
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
    parser.add_argument("--c2-initial", type=float, default=0.0)
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

    summary = []
    traces = {}
    for irradiance in WANG_INTENSITY_VALUES_MW_MM2:
        trace = simulate_diagnostic(args.cell_params, rows, gbar, irradiance, args, protocol)
        traces[irradiance] = trace
        peak, _, ttp, tau = current_metrics(trace["t"], trace["clamp_i"], protocol)
        summary.append(
            {
                "label": args.label,
                "irradiance_mw_mm2": irradiance,
                "gbar_mS_cm2": gbar,
                "clamp_peak_nA": peak,
                "clamp_time_to_peak_ms": ttp,
                "clamp_tau_ms": tau,
                "local_peak_nA": float(np.max(trace["local_i"])),
                "local_time_to_peak_ms": first_peak_time(trace["t"], trace["local_i"]),
                "open_peak": float(np.max(trace["open"])),
                "open_time_to_peak_ms": first_peak_time(trace["t"], trace["open"]),
                "o1_peak": float(np.max(trace["o1"])),
                "o1_time_to_peak_ms": first_peak_time(trace["t"], trace["o1"]),
                "o2_peak": float(np.max(trace["o2"])),
                "o2_time_to_peak_ms": first_peak_time(trace["t"], trace["o2"]),
                "c2_at_100ms": float(trace["c2"][np.argmin(np.abs(trace["t"] - 100.0))]),
                "p_time_to_90pct_ms": float(trace["t"][np.argmax(trace["p"] >= 0.9)]) if np.any(trace["p"] >= 0.9) else float("nan"),
            }
        )

    with (args.out_dir / "intensity_state_diagnostics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(WANG_INTENSITY_VALUES_MW_MM2)))
    for color, irradiance in zip(colors, WANG_INTENSITY_VALUES_MW_MM2):
        trace = traces[irradiance]
        label = f"{irradiance:g}"
        axes[0, 0].plot(trace["t"], trace["clamp_i"], color=color, lw=1.1, label=label)
        axes[0, 1].plot(trace["t"], trace["open"], color=color, lw=1.1)
        axes[1, 0].plot(trace["t"], trace["o1"], color=color, lw=1.1)
        axes[1, 1].plot(trace["t"], trace["o2"], color=color, lw=1.1)
    axes[0, 0].set_title("SEClamp current")
    axes[0, 1].set_title("weighted open = O1 + gamma O2")
    axes[1, 0].set_title("weighted O1")
    axes[1, 1].set_title("weighted O2")
    for ax in axes.flat:
        ax.set_xlim(-2, 80)
        ax.set_xlabel("time (ms)")
    axes[0, 0].set_ylabel("nA")
    axes[0, 1].set_ylabel("fraction")
    axes[1, 0].set_ylabel("fraction")
    axes[1, 1].set_ylabel("fraction")
    axes[0, 0].legend(title="mW/mm2", fontsize=8)
    fig.savefig(args.out_dir / "intensity_state_traces.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    xs = [row["irradiance_mw_mm2"] for row in summary]
    ax.plot(xs, [row["clamp_time_to_peak_ms"] for row in summary], marker="o", label="SEClamp current")
    ax.plot(xs, [row["local_time_to_peak_ms"] for row in summary], marker="o", label="summed local ChR2 current")
    ax.plot(xs, [row["open_time_to_peak_ms"] for row in summary], marker="o", label="weighted open state")
    ax.plot(xs, [row["o1_time_to_peak_ms"] for row in summary], marker="o", label="weighted O1")
    ax.set_xscale("log")
    ax.set_xlabel("irradiance (mW/mm2)")
    ax.set_ylabel("time to peak (ms)")
    ax.legend()
    fig.savefig(args.out_dir / "time_to_peak_by_intensity.png", dpi=220)
    plt.close(fig)

    for row in summary:
        print(row)


if __name__ == "__main__":
    main()
