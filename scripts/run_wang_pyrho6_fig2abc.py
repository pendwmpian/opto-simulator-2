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
from run_wang_electrode_amplifier_model import apply_electrode_amplifier_readout
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    add_localization_weights,
    fit_hill_k,
    load_mechanisms,
    prepare_template_rows,
    soma_section,
)


def parse_floats(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def q10_rate_scale(q10: float, reference_c: float, target_c: float) -> float:
    if q10 <= 0:
        raise ValueError("Q10 must be positive")
    return 1.0 / (q10 ** ((reference_c - target_c) / 10.0))


def install_chr2_pyrho6(
    cell: CellFromNetPyNE,
    rows_template,
    gbar_mS_cm2: float,
    irradiance_scale: float,
    rate_scale: float,
    gamma_scale: float,
    phi_m_scale: float,
):
    for sec in cell.sections.values():
        sec.insert("chr2_pyrho6")
        for seg in sec:
            seg.chr2_pyrho6.gbar = 0.0
            seg.chr2_pyrho6.irr = 0.0
            seg.chr2_pyrho6.q10_scale = rate_scale
            seg.chr2_pyrho6.gamma *= gamma_scale
            seg.chr2_pyrho6.phi_m_scale = phi_m_scale

    rows = []
    for template in rows_template:
        sec_name = template["sec"].name().split(".")[-1]
        sec = cell.sections[sec_name]
        seg = sec(template["xloc"])
        seg.chr2_pyrho6.gbar = gbar_mS_cm2 * template.get("localization_weight", 1.0)
        seg.chr2_pyrho6.q10_scale = rate_scale
        seg.chr2_pyrho6.gamma *= gamma_scale
        seg.chr2_pyrho6.phi_m_scale = phi_m_scale
        rows.append(
            {
                "seg": seg,
                "base_irradiance_mw_mm2": template["irradiance_mw_mm2"] * irradiance_scale,
                "area_um2": template["area_um2"],
            }
        )
    return rows


def simulate_pyrho6_current(
    cell_params: Path,
    rows_template,
    gbar_mS_cm2: float,
    irradiance_scale: float,
    rs_mohm: float,
    protocol: Wang2007Protocol,
    dt_ms: float,
    tstop_ms: float,
    pulse_duration_ms: float,
    protocol_irradiance_scale: float,
    rate_scale: float,
    gamma_scale: float,
    phi_m_scale: float,
    use_original_biophysics: bool,
    subtract_no_light: bool = True,
    return_raw_current: bool = False,
):
    if subtract_no_light:
        t_light, i_light, v_light = simulate_pyrho6_current(
            cell_params,
            rows_template,
            gbar_mS_cm2,
            irradiance_scale,
            rs_mohm,
            protocol,
            dt_ms,
            tstop_ms,
            pulse_duration_ms,
            protocol_irradiance_scale,
            rate_scale,
            gamma_scale,
            phi_m_scale,
            use_original_biophysics,
            subtract_no_light=False,
            return_raw_current=True,
        )
        t_base, i_base, _ = simulate_pyrho6_current(
            cell_params,
            rows_template,
            gbar_mS_cm2,
            irradiance_scale,
            rs_mohm,
            protocol,
            dt_ms,
            tstop_ms,
            0.0,
            0.0,
            rate_scale,
            gamma_scale,
            phi_m_scale,
            use_original_biophysics,
            subtract_no_light=False,
            return_raw_current=True,
        )
        if len(i_base) != len(i_light):
            i_base = np.interp(t_light, t_base, i_base)
        return t_light, np.abs(i_light - i_base), v_light

    h.load_file("stdrun.hoc")
    h.celsius = 22.0
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=use_original_biophysics)
    rows = install_chr2_pyrho6(cell, rows_template, gbar_mS_cm2, irradiance_scale, rate_scale, gamma_scale, phi_m_scale)

    clamp = h.SEClamp(soma_section(cell)(0.5))
    clamp.dur1 = tstop_ms + 1.0
    clamp.amp1 = protocol.holding_voltage_mV
    clamp.rs = rs_mohm

    h.dt = dt_ms
    h.tstop = tstop_ms
    h.finitialize(protocol.holding_voltage_mV)
    times = []
    currents = []
    voltages = []
    soma = soma_section(cell)(0.5)
    while h.t <= tstop_ms:
        light_on = 0.0 <= h.t < pulse_duration_ms
        for row in rows:
            row["seg"].chr2_pyrho6.irr = row["base_irradiance_mw_mm2"] * protocol_irradiance_scale if light_on else 0.0
        times.append(float(h.t))
        currents.append(float(clamp.i))
        voltages.append(float(soma.v))
        h.fadvance()
    current = np.asarray(currents)
    if not return_raw_current:
        current = np.abs(current)
    return np.asarray(times), current, np.asarray(voltages)


def calibrate_pyrho6_gbar(cell_params, rows, args, protocol, rate_scale: float, gamma_scale: float, phi_m_scale: float):
    lo = 1.0e-5
    hi = 2.0
    last = None
    for _ in range(args.binary_iterations):
        mid = (lo * hi) ** 0.5
        t, current, _ = simulate_pyrho6_current(
            cell_params,
            rows,
            mid,
            args.irradiance_scale,
            args.rs_mohm,
            protocol,
            args.dt_ms,
            max(180.0, args.saturating_duration_ms + 80.0),
            args.saturating_duration_ms,
            args.saturating_irradiance_mw_mm2 / protocol.irradiance_mw_mm2,
            rate_scale,
            gamma_scale,
            phi_m_scale,
            args.use_original_biophysics,
            subtract_no_light=args.use_original_biophysics,
        )
        peak = float(np.max(current))
        last = (mid, peak)
        if peak < protocol.target_imax_nA:
            lo = mid
        else:
            hi = mid
    if last is None:
        raise RuntimeError("gbar calibration did not run")
    return last


def run_trace(cell_params, rows, gbar, irradiance, pulse_duration_ms, tstop_ms, args, protocol, rate_scale, gamma_scale, phi_m_scale):
    t, clamp_i, soma_v = simulate_pyrho6_current(
        cell_params,
        rows,
        gbar,
        args.irradiance_scale,
        args.rs_mohm,
        protocol,
        args.dt_ms,
        tstop_ms,
        pulse_duration_ms,
        irradiance / protocol.irradiance_mw_mm2,
        rate_scale,
        gamma_scale,
        phi_m_scale,
        args.use_original_biophysics,
        subtract_no_light=args.use_original_biophysics,
    )
    filtered_i = apply_electrode_amplifier_readout(t, clamp_i, args)
    return {"t": t, "clamp_i": clamp_i, "filtered_i": filtered_i, "soma_v": soma_v}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang_pyrho6_fig2abc"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--q10-values", default="1.0,1.5,2.0")
    parser.add_argument("--gamma-scales", default="1.0")
    parser.add_argument("--phi-m-scales", default="1.0")
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

    summary_rows = []
    intensity_rows = []
    traces = {}
    for q10 in parse_floats(args.q10_values):
        rate_scale = q10_rate_scale(q10, args.reference_c, args.target_c)
        for gamma_scale in parse_floats(args.gamma_scales):
            for phi_m_scale in parse_floats(args.phi_m_scales):
                gbar, calib_peak = calibrate_pyrho6_gbar(args.cell_params, rows, args, protocol, rate_scale, gamma_scale, phi_m_scale)
                trace_1s = run_trace(args.cell_params, rows, gbar, protocol.irradiance_mw_mm2, 1000.0, 1200.0, args, protocol, rate_scale, gamma_scale, phi_m_scale)
                traces[(q10, gamma_scale, phi_m_scale)] = trace_1s
                raw_peak, raw_steady, raw_ttp, raw_tau = current_metrics(trace_1s["t"], trace_1s["clamp_i"], protocol)
                filt_peak, filt_steady, filt_ttp, filt_tau = current_metrics(trace_1s["t"], trace_1s["filtered_i"], protocol)

                raw_peaks = []
                filt_peaks = []
                for irr in WANG_INTENSITY_VALUES_MW_MM2:
                    trace = run_trace(args.cell_params, rows, gbar, irr, 100.0, 180.0, args, protocol, rate_scale, gamma_scale, phi_m_scale)
                    raw_peak_i, raw_steady_i, raw_ttp_i, raw_tau_i = current_metrics(trace["t"], trace["clamp_i"], protocol)
                    filt_peak_i, filt_steady_i, filt_ttp_i, filt_tau_i = current_metrics(trace["t"], trace["filtered_i"], protocol)
                    raw_peaks.append(raw_peak_i)
                    filt_peaks.append(filt_peak_i)
                    intensity_rows.append(
                        {
                            "q10": q10,
                            "gamma_scale": gamma_scale,
                            "phi_m_scale": phi_m_scale,
                            "effective_gamma": 0.00369 * gamma_scale,
                            "effective_phi_m": 5.02e17 * phi_m_scale,
                            "rate_scale": rate_scale,
                            "irradiance_mw_mm2": irr,
                            "gbar_mS_cm2": gbar,
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
                filt_k, filt_imax, filt_hill = fit_hill_k(WANG_INTENSITY_VALUES_MW_MM2, filt_peaks)
                summary_rows.append(
                    {
                        "label": args.label,
                        "model": "PyRhO6",
                        "q10": q10,
                        "gamma_scale": gamma_scale,
                        "phi_m_scale": phi_m_scale,
                        "effective_gamma": 0.00369 * gamma_scale,
                        "effective_phi_m": 5.02e17 * phi_m_scale,
                        "rate_scale": rate_scale,
                        "gbar_mS_cm2": gbar,
                        "calibration_peak_nA": calib_peak,
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
                    }
                )

    with (args.out_dir / "pyrho6_fig2abc_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    with (args.out_dir / "pyrho6_fig2abc_intensity.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(intensity_rows[0].keys()))
        writer.writeheader()
        writer.writerows(intensity_rows)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2), constrained_layout=True)
    for (q10, gamma_scale, phi_m_scale), trace in traces.items():
        axes[0].plot(trace["t"], -trace["filtered_i"], lw=1.2, label=f"Q10={q10:g}, gamma x{gamma_scale:g}, phi x{phi_m_scale:g}")
    axes[0].set_xlim(-10, 1100)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("current (nA)")
    axes[0].set_title("Fig. 2A-like 9.2 mW/mm2")
    axes[0].legend(fontsize=8)

    for row in summary_rows:
        q10 = row["q10"]
        gamma_scale = row["gamma_scale"]
        phi_m_scale = row["phi_m_scale"]
        sub = [
            item
            for item in intensity_rows
            if item["q10"] == q10 and item["gamma_scale"] == gamma_scale and item["phi_m_scale"] == phi_m_scale
        ]
        axes[1].plot(
            [item["irradiance_mw_mm2"] for item in sub],
            [item["filtered_peak_nA"] for item in sub],
            marker="o",
            label=f"Q10={q10:g}, g x{gamma_scale:g}, phi x{phi_m_scale:g}, K={row['filtered_intensity_k_mw_mm2']:.2g}",
        )
    axes[1].axvline(0.84, color="#d62728", ls="--", lw=1, label="Wang K~0.84")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("peak current (nA)")
    axes[1].set_title("Fig. 2B/C-like intensity response")
    axes[1].legend(fontsize=7)

    x = np.arange(len(summary_rows))
    axes[2].bar(x - 0.18, [row["filtered_peak_9p2_1s_nA"] for row in summary_rows], width=0.36, label="peak")
    axes[2].bar(x + 0.18, [row["filtered_steady_9p2_1s_nA"] for row in summary_rows], width=0.36, label="1 s steady")
    axes[2].axhline(protocol.target_peak_current_nA, color="#777777", ls=":", lw=1, label="Wang 9.2 peak")
    axes[2].axhline(0.4, color="#d62728", ls="--", lw=1, label="Wang visual steady ~0.4")
    axes[2].set_xticks(x, [f"Q{row['q10']:g}\ng x{row['gamma_scale']:g}\nphi x{row['phi_m_scale']:g}" for row in summary_rows])
    axes[2].set_ylabel("current (nA)")
    axes[2].set_title("Peak and sustained current")
    axes[2].legend(fontsize=7)
    fig.savefig(args.out_dir / "pyrho6_fig2abc_summary.png", dpi=220)
    plt.close(fig)

    for row in summary_rows:
        print(row)


if __name__ == "__main__":
    main()
