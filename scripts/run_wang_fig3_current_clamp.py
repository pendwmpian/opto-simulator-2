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
from scipy.optimize import curve_fit

from calibrate_wang2007_gbar import Wang2007Protocol
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_electrode_amplifier_model import first_order_lowpass
from run_wang_seclamp_mod_suite import (
    WANG_INTENSITY_VALUES_MW_MM2,
    add_localization_weights,
    hill_response,
    install_chr2,
    load_mechanisms,
    prepare_template_rows,
    soma_section,
    williams_q10_scales,
)


def spike_peaks(t: np.ndarray, v: np.ndarray, threshold_mV: float = 0.0, refractory_ms: float = 2.0):
    peaks = []
    last_t = -1.0e9
    for i in range(1, len(v) - 1):
        if t[i] < 0:
            continue
        if v[i] >= threshold_mV and v[i] >= v[i - 1] and v[i] > v[i + 1] and t[i] - last_t >= refractory_ms:
            peaks.append((float(t[i]), float(v[i])))
            last_t = float(t[i])
    return peaks


def fit_spike_count_k(intensities: list[float], counts: list[float]):
    xs = np.asarray(intensities, dtype=float)
    ys = np.asarray(counts, dtype=float)
    if np.max(ys) <= 0:
        return float("nan"), float("nan"), float("nan")
    try:
        popt, _ = curve_fit(
            hill_response,
            xs,
            ys,
            p0=[float(np.max(ys)), 0.49, 0.82],
            bounds=([1.0e-9, 1.0e-6, 0.1], [np.inf, np.inf, 5.0]),
            maxfev=20000,
        )
    except RuntimeError:
        return float("nan"), float(np.max(ys)), float("nan")
    imax, k_half, hill_n = [float(value) for value in popt]
    return k_half, imax, hill_n


def make_pulse_train(duration_ms: float, frequency_hz: float, train_duration_ms: float):
    interval_ms = 1000.0 / frequency_hz
    starts = np.arange(0.0, train_duration_ms, interval_ms)
    return [(float(start), float(start + duration_ms)) for start in starts]


def light_is_on(t_ms: float, pulses: list[tuple[float, float]]) -> bool:
    return any(start <= t_ms < stop for start, stop in pulses)


def simulate_current_clamp(
    cell_params: Path,
    rows_template,
    gbar_mS_cm2: float,
    irradiance_mw_mm2: float,
    pulses: list[tuple[float, float]],
    tstop_ms: float,
    args,
    protocol: Wang2007Protocol,
):
    h.load_file("stdrun.hoc")
    h.celsius = args.target_c
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=args.use_original_biophysics)
    rows = install_chr2(
        cell,
        rows_template,
        gbar_mS_cm2,
        args.irradiance_scale,
        williams_q10_scales(args.reference_c, args.target_c),
        True,
    )

    stim = h.IClamp(soma_section(cell)(0.5))
    stim.delay = -1.0e9
    stim.dur = 1.0e9
    stim.amp = args.dc_bias_nA

    soma = soma_section(cell)(0.5)
    h.dt = args.dt_ms
    h.tstop = tstop_ms
    h.finitialize(args.initial_v_mV)

    times = []
    voltages = []
    light_flags = []
    start_t = -args.equilibration_ms
    h.t = start_t
    while h.t <= tstop_ms:
        on = light_is_on(float(h.t), pulses)
        scale = irradiance_mw_mm2 / protocol.irradiance_mw_mm2
        for row in rows:
            row["seg"].chr2_4state.irr = row["base_irradiance_mw_mm2"] * scale if on else 0.0
        if h.t >= 0.0:
            times.append(float(h.t))
            voltages.append(float(soma.v))
            light_flags.append(1 if on else 0)
        h.fadvance()

    t = np.asarray(times)
    v = np.asarray(voltages)
    displayed_v = first_order_lowpass(t, v, args.voltage_filter_tau_ms) if args.voltage_filter_tau_ms > 0 else v.copy()
    return {
        "t": t,
        "v": v,
        "displayed_v": displayed_v,
        "light": np.asarray(light_flags),
        "rest_v_mV": float(v[0]) if len(v) else float("nan"),
    }


def run_intensity_series(rows, args, protocol):
    summaries = []
    traces = {}
    for irr in WANG_INTENSITY_VALUES_MW_MM2:
        trace = simulate_current_clamp(
            args.cell_params,
            rows,
            args.gbar_mS_cm2,
            irr,
            [(0.0, 1000.0)],
            1100.0,
            args,
            protocol,
        )
        spikes = spike_peaks(trace["t"], trace["v"], args.spike_threshold_mV, args.refractory_ms)
        first_latency = spikes[0][0] if spikes else float("nan")
        summaries.append(
            {
                "irradiance_mw_mm2": irr,
                "rest_v_mV": trace["rest_v_mV"],
                "spike_count_1s": len(spikes),
                "first_ap_peak_latency_ms": first_latency,
                "max_v_mV": float(np.max(trace["v"])),
            }
        )
        traces[irr] = trace
    k_half, imax, hill_n = fit_spike_count_k(
        [row["irradiance_mw_mm2"] for row in summaries],
        [row["spike_count_1s"] for row in summaries],
    )
    return summaries, traces, {"spike_count_k_mw_mm2": k_half, "spike_count_imax": imax, "spike_count_hill_n": hill_n}


def run_frequency_following(rows, args, protocol):
    summaries = []
    traces = {}
    for freq in args.frequency_values_hz:
        pulses = make_pulse_train(4.0, freq, args.frequency_train_duration_ms)
        trace = simulate_current_clamp(
            args.cell_params,
            rows,
            args.gbar_mS_cm2,
            9.2,
            pulses,
            args.frequency_train_duration_ms + 100.0,
            args,
            protocol,
        )
        spikes = spike_peaks(trace["t"], trace["v"], args.spike_threshold_mV, args.refractory_ms)
        successes = 0
        latencies = []
        for start, _ in pulses:
            candidates = [spike_t for spike_t, _ in spikes if start <= spike_t < start + args.evoked_window_ms]
            if candidates:
                successes += 1
                latencies.append(candidates[0] - start)
        summaries.append(
            {
                "frequency_hz": freq,
                "pulse_count": len(pulses),
                "evoked_spike_count": successes,
                "success_probability": successes / len(pulses) if pulses else float("nan"),
                "mean_latency_ms": float(np.mean(latencies)) if latencies else float("nan"),
                "rest_v_mV": trace["rest_v_mV"],
            }
        )
        traces[freq] = trace
    return summaries, traces


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_intensity(out_dir: Path, summaries, traces, fit):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    for irr in [0.58, 1.15, 2.3, 9.2]:
        if irr in traces:
            trace = traces[irr]
            axes[0].plot(trace["t"], trace["displayed_v"], lw=1.0, label=f"{irr:g}")
    axes[0].set_xlim(-10, 220)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("soma Vm (mV)")
    axes[0].set_title("1 s flash examples")
    axes[0].legend(title="mW/mm2", fontsize=8)

    xs = [row["irradiance_mw_mm2"] for row in summaries]
    counts = [row["spike_count_1s"] for row in summaries]
    axes[1].plot(xs, counts, marker="o")
    axes[1].axvline(0.49, color="#d62728", ls="--", lw=1)
    if any(x > 0 for x in xs):
        axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("spikes / 1 s")
    axes[1].set_title(f"K={fit['spike_count_k_mw_mm2']:.3g}")

    latencies = [row["first_ap_peak_latency_ms"] for row in summaries]
    if any(np.isfinite(latencies)):
        axes[2].plot(xs, latencies, marker="o")
    else:
        axes[2].text(0.5, 0.5, "no spikes", transform=axes[2].transAxes, ha="center", va="center")
    axes[2].axhline(6.1, color="#d62728", ls="--", lw=1)
    if any(x > 0 for x in xs):
        axes[2].set_xscale("log")
    axes[2].set_xlabel("irradiance (mW/mm2)")
    axes[2].set_ylabel("first AP peak latency (ms)")
    axes[2].set_title("first spike latency")
    fig.savefig(out_dir / "fig3_intensity_spike_count.png", dpi=220)
    plt.close(fig)


def plot_frequency(out_dir: Path, summaries, traces):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), constrained_layout=True)
    for freq in [10.0, 30.0, 40.0]:
        if freq in traces:
            trace = traces[freq]
            axes[0].plot(trace["t"], trace["displayed_v"], lw=1.0, label=f"{freq:g} Hz")
    axes[0].set_xlim(-10, 260)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("soma Vm (mV)")
    axes[0].set_title("4 ms pulse trains")
    axes[0].legend(fontsize=8)

    xs = [row["frequency_hz"] for row in summaries]
    probs = [row["success_probability"] for row in summaries]
    axes[1].plot(xs, probs, marker="o")
    axes[1].axvline(30.0, color="#d62728", ls="--", lw=1)
    axes[1].axhline(0.5, color="#999999", ls=":", lw=1)
    axes[1].set_xlabel("frequency (Hz)")
    axes[1].set_ylabel("evoked spike probability")
    axes[1].set_ylim(-0.05, 1.05)
    fig.savefig(out_dir / "fig3_frequency_following.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_fig3_current_clamp"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--gbar-mS-cm2", type=float, default=0.13478410273935307)
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
    parser.add_argument("--initial-v-mV", type=float, default=-70.0)
    parser.add_argument("--equilibration-ms", type=float, default=500.0)
    parser.add_argument("--dc-bias-nA", type=float, default=0.0)
    parser.add_argument("--dt-ms", type=float, default=0.05)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    parser.add_argument("--refractory-ms", type=float, default=2.0)
    parser.add_argument("--voltage-filter-tau-ms", type=float, default=0.0)
    parser.add_argument("--use-original-biophysics", action="store_true")
    parser.add_argument("--frequency-values-hz", type=lambda s: [float(x) for x in s.split(",")], default="5,10,20,30,35,40,50")
    parser.add_argument("--frequency-train-duration-ms", type=float, default=1000.0)
    parser.add_argument("--evoked-window-ms", type=float, default=20.0)
    parser.add_argument("--skip-intensity", action="store_true")
    parser.add_argument("--skip-frequency", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    protocol = Wang2007Protocol()
    h.load_file("stdrun.hoc")
    template_cell = CellFromNetPyNE(args.cell_params, use_original_biophysics=args.use_original_biophysics)
    rows = prepare_template_rows(template_cell, protocol, args)
    rows = add_localization_weights(rows, args.localization, args.proximal_cutoff_um)

    if not args.skip_intensity:
        intensity_rows, intensity_traces, fit = run_intensity_series(rows, args, protocol)
        write_csv(args.out_dir / "fig3_intensity_spike_count.csv", intensity_rows)
        with (args.out_dir / "fig3_intensity_fit.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(fit.keys()))
            writer.writeheader()
            writer.writerow(fit)
        plot_intensity(args.out_dir, intensity_rows, intensity_traces, fit)
        print(fit)
        for row in intensity_rows:
            print(row)
    if not args.skip_frequency:
        frequency_rows, frequency_traces = run_frequency_following(rows, args, protocol)
        write_csv(args.out_dir / "fig3_frequency_following.csv", frequency_rows)
        plot_frequency(args.out_dir, frequency_rows, frequency_traces)
        for row in frequency_rows:
            print(row)


if __name__ == "__main__":
    main()
