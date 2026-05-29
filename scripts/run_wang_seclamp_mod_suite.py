from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from neuron import h

from calibrate_wang2007_gbar import (
    Wang2007Protocol,
    add_vitro_irradiance,
    current_metrics,
    half_max_x,
    segment_rows_vertical_slice,
    soma_axis_value,
)
from run_chr2_single_cell_sweep import CellFromNetPyNE

WANG_INTENSITY_VALUES_MW_MM2 = [0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2]


@dataclass(frozen=True)
class ModSuiteResult:
    label: str
    localization: str
    calibration_target: str
    calibration_irradiance_mw_mm2: float
    calibration_duration_ms: float
    calibration_target_current_nA: float
    rs_mohm: float
    irradiance_scale: float
    q10: float
    q10_rate_scale: float
    q10_mode: str
    gbar_mS_cm2: float
    peak_1s_9p2_nA: float
    time_to_peak_1s_9p2_ms: float
    inactivation_tau_1s_9p2_ms: float | None
    soma_v_min_mV: float
    soma_v_max_mV: float
    soma_v_mean_mV: float
    intensity_k_mw_mm2: float | None
    intensity_imax_nA: float
    intensity_hill_n: float
    duration_k_ms: float | None
    duration_imax_nA: float
    retained_area_um2: float
    weighted_area_um2: float
    mean_irradiance_mw_mm2: float


def hill_response(x, imax, k_half, hill_n):
    x = np.asarray(x, dtype=float)
    return imax * (x**hill_n) / (k_half**hill_n + x**hill_n)


def fit_hill_k(xs: list[float], ys: list[float]) -> tuple[float | None, float, float]:
    xs_arr = np.asarray(xs, dtype=float)
    ys_arr = np.asarray(ys, dtype=float)
    if len(xs_arr) < 4 or float(np.max(ys_arr)) <= 0.0:
        return None, float(np.max(ys_arr)), float("nan")
    try:
        popt, _ = curve_fit(
            hill_response,
            xs_arr,
            ys_arr,
            p0=[float(np.max(ys_arr)), 0.84, 0.76],
            bounds=([1.0e-9, 1.0e-6, 0.1], [np.inf, np.inf, 5.0]),
            maxfev=20000,
        )
    except RuntimeError:
        return half_max_x(xs, ys), float(np.max(ys_arr)), float("nan")
    imax, k_half, hill_n = [float(value) for value in popt]
    return k_half, imax, hill_n


def load_mechanisms(repo_root: Path) -> None:
    lib = repo_root / "arm64" / ".libs" / "libnrnmech.dylib"
    if not lib.exists():
        lib = repo_root / "arm64" / "libnrnmech.dylib"
    try:
        h.nrn_load_dll(str(lib))
    except RuntimeError as exc:
        if "already exists" not in str(exc):
            raise


def soma_section(cell: CellFromNetPyNE):
    soma_name = "soma" if "soma" in cell.sections else next(name for name in cell.sections if name.startswith("soma"))
    return cell.sections[soma_name]


def prepare_template_rows(cell: CellFromNetPyNE, protocol: Wang2007Protocol, args):
    soma_axis = soma_axis_value(cell, args.slice_axis)
    slice_center_um = soma_axis + protocol.slice_thickness_um / 2.0 - args.soma_entry_depth_um
    illumination_entry_um = soma_axis - args.soma_entry_depth_um
    rows, _ = segment_rows_vertical_slice(
        cell,
        args.slice_axis,
        slice_center_um,
        protocol.slice_thickness_um,
        args.truncate_radius_um,
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


def add_localization_weights(rows, mode: str, proximal_cutoff_um: float = 150.0):
    soma_rows = [row for row in rows if row["kind"] == "soma"]
    if soma_rows:
        soma_xyz = np.array([soma_rows[0]["x"], soma_rows[0]["y"], soma_rows[0]["z"]], dtype=float)
    else:
        soma_xyz = np.array([rows[0]["x"], rows[0]["y"], rows[0]["z"]], dtype=float)
    weighted = []
    for row in rows:
        kind = row["kind"]
        xyz = np.array([row["x"], row["y"], row["z"]], dtype=float)
        dist = float(np.linalg.norm(xyz - soma_xyz))
        is_soma = kind == "soma"
        is_dend = kind in {"apic", "dend"}
        is_prox = is_dend and dist <= proximal_cutoff_um
        is_dist = is_dend and dist > proximal_cutoff_um
        if mode == "uniform":
            weight = 1.0
        elif mode == "soma-only":
            weight = 1.0 if is_soma else 0.0
        elif mode == "dendrite-only":
            weight = 1.0 if is_dend else 0.0
        elif mode == "proximal-enriched":
            weight = 1.0 if (is_soma or is_prox) else (0.25 if is_dist else 0.0)
        elif mode == "proximal-reduced":
            weight = 0.25 if (is_soma or is_prox) else (1.0 if is_dist else 0.0)
        else:
            raise ValueError(f"Unknown localization mode: {mode}")
        weighted.append({**row, "localization_weight": weight, "distance_from_soma_um": dist})
    return weighted


def williams_q10_scales(reference_c: float = 37.0, target_c: float = 22.0) -> dict[str, float]:
    exponent = (reference_c - target_c) / 10.0
    q10 = {
        "k1": 1.46,
        "k2": 2.77,
        "gd1": 1.97,
        "gd2": 1.77,
        "e12": 1.10,
        "e21": 1.95,
        "gr": 2.56,
    }
    return {name: 1.0 / (value**exponent) for name, value in q10.items()}


def generic_q10_scales(q10: float, reference_c: float = 37.0, target_c: float = 22.0) -> dict[str, float]:
    scale = 1.0 / (q10 ** ((reference_c - target_c) / 10.0))
    return {name: scale for name in ["k1", "k2", "gd1", "gd2", "e12", "e21", "gr"]}


def install_chr2(
    cell: CellFromNetPyNE,
    rows_template,
    gbar_mS_cm2: float,
    irradiance_scale: float,
    q10_scales: dict[str, float],
    use_photon_flux: bool,
):
    for sec in cell.sections.values():
        sec.insert("chr2_4state")
        for seg in sec:
            seg.chr2_4state.gbar = 0.0
            seg.chr2_4state.irr = 0.0
            seg.chr2_4state.activation_mode = 1.0 if use_photon_flux else 0.0
            seg.chr2_4state.q10_scale = 1.0
            seg.chr2_4state.q10_k1_scale = q10_scales["k1"]
            seg.chr2_4state.q10_k2_scale = q10_scales["k2"]
            seg.chr2_4state.q10_gd1_scale = q10_scales["gd1"]
            seg.chr2_4state.q10_gd2_scale = q10_scales["gd2"]
            seg.chr2_4state.q10_e12_scale = q10_scales["e12"]
            seg.chr2_4state.q10_e21_scale = q10_scales["e21"]
            seg.chr2_4state.q10_gr_scale = q10_scales["gr"]

    rows = []
    for template in rows_template:
        sec_name = template["sec"].name().split(".")[-1]
        sec = cell.sections[sec_name]
        seg = sec(template["xloc"])
        seg.chr2_4state.gbar = gbar_mS_cm2 * template.get("localization_weight", 1.0)
        seg.chr2_4state.activation_mode = 1.0 if use_photon_flux else 0.0
        seg.chr2_4state.q10_scale = 1.0
        seg.chr2_4state.q10_k1_scale = q10_scales["k1"]
        seg.chr2_4state.q10_k2_scale = q10_scales["k2"]
        seg.chr2_4state.q10_gd1_scale = q10_scales["gd1"]
        seg.chr2_4state.q10_gd2_scale = q10_scales["gd2"]
        seg.chr2_4state.q10_e12_scale = q10_scales["e12"]
        seg.chr2_4state.q10_e21_scale = q10_scales["e21"]
        seg.chr2_4state.q10_gr_scale = q10_scales["gr"]
        rows.append(
            {
                "seg": seg,
                "base_irradiance_mw_mm2": template["irradiance_mw_mm2"] * irradiance_scale,
                "area_um2": template["area_um2"],
            }
        )
    return rows


def simulate_mod_current(
    cell_params: Path,
    rows_template,
    gbar_mS_cm2: float,
    irradiance_scale: float,
    rs_mohm: float,
    protocol: Wang2007Protocol,
    dt_ms: float,
    tstop_ms: float,
    pulse_duration_ms: float,
    protocol_irradiance_scale: float = 1.0,
    q10_scales: dict[str, float] | None = None,
    use_photon_flux: bool = False,
    use_original_biophysics: bool = False,
    subtract_no_light: bool = False,
    return_raw_current: bool = False,
):
    if subtract_no_light:
        t_light, i_light, v_light = simulate_mod_current(
            cell_params,
            rows_template,
            gbar_mS_cm2,
            irradiance_scale,
            rs_mohm,
            protocol,
            dt_ms,
            tstop_ms,
            pulse_duration_ms,
            protocol_irradiance_scale=protocol_irradiance_scale,
            q10_scales=q10_scales,
            use_photon_flux=use_photon_flux,
            use_original_biophysics=use_original_biophysics,
            subtract_no_light=False,
            return_raw_current=True,
        )
        t_base, i_base, _ = simulate_mod_current(
            cell_params,
            rows_template,
            gbar_mS_cm2,
            irradiance_scale,
            rs_mohm,
            protocol,
            dt_ms,
            tstop_ms,
            0.0,
            protocol_irradiance_scale=0.0,
            q10_scales=q10_scales,
            use_photon_flux=use_photon_flux,
            use_original_biophysics=use_original_biophysics,
            subtract_no_light=False,
            return_raw_current=True,
        )
        if len(i_base) != len(i_light):
            i_base = np.interp(t_light, t_base, i_base)
        return t_light, np.abs(i_light - i_base), v_light

    h.load_file("stdrun.hoc")
    h.celsius = 22.0
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=use_original_biophysics)
    rows = install_chr2(
        cell,
        rows_template,
        gbar_mS_cm2,
        irradiance_scale,
        q10_scales or generic_q10_scales(1.0),
        use_photon_flux,
    )

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
            row["seg"].chr2_4state.irr = row["base_irradiance_mw_mm2"] * protocol_irradiance_scale if light_on else 0.0
        times.append(float(h.t))
        currents.append(float(clamp.i))
        voltages.append(float(soma.v))
        h.fadvance()
    current = np.asarray(currents)
    if not return_raw_current:
        current = np.abs(current)
    return np.asarray(times), current, np.asarray(voltages)


def calibrate_gbar(
    cell_params,
    rows,
    irradiance_scale,
    rs_mohm,
    protocol,
    dt_ms,
    iterations,
    q10_scales,
    use_photon_flux,
    target_current_nA,
    calibration_irradiance_mw_mm2,
    calibration_duration_ms,
    use_original_biophysics=False,
):
    lo = 1.0e-5
    hi = 1.0
    last = None
    for _ in range(iterations):
        mid = (lo * hi) ** 0.5
        t, current, _ = simulate_mod_current(
            cell_params,
            rows,
            mid,
            irradiance_scale,
            rs_mohm,
            protocol,
            dt_ms,
            max(180.0, calibration_duration_ms + 80.0),
            calibration_duration_ms,
            protocol_irradiance_scale=calibration_irradiance_mw_mm2 / protocol.irradiance_mw_mm2,
            q10_scales=q10_scales,
            use_photon_flux=use_photon_flux,
            use_original_biophysics=use_original_biophysics,
            subtract_no_light=use_original_biophysics,
        )
        peak = float(np.max(current))
        last = (mid, t, current)
        if peak < target_current_nA:
            lo = mid
        else:
            hi = mid
    if last is None:
        raise RuntimeError("Calibration did not run")
    return last


def run_condition(label, cell_params, rows, rs_mohm, irradiance_scale, q10, protocol, args):
    rows = add_localization_weights(rows, args.localization, args.proximal_cutoff_um)
    if args.q10_mode == "williams":
        q10_scales = williams_q10_scales(args.reference_c, args.target_c)
        q10_rate_scale = float("nan")
    else:
        q10_scales = generic_q10_scales(q10, args.reference_c, args.target_c)
        q10_rate_scale = q10_scales["k1"]
    use_photon_flux = args.activation_mode == "photon"
    if args.calibration_target == "imax_saturating":
        target_current_nA = protocol.target_imax_nA
        calibration_irradiance_mw_mm2 = args.saturating_irradiance_mw_mm2
        calibration_duration_ms = args.saturating_duration_ms
    else:
        target_current_nA = protocol.target_peak_current_nA
        calibration_irradiance_mw_mm2 = protocol.irradiance_mw_mm2
        calibration_duration_ms = protocol.pulse_duration_ms

    gbar, _, _ = calibrate_gbar(
        cell_params,
        rows,
        irradiance_scale,
        rs_mohm,
        protocol,
        args.dt_ms,
        args.binary_iterations,
        q10_scales,
        use_photon_flux,
        target_current_nA,
        calibration_irradiance_mw_mm2,
        calibration_duration_ms,
        getattr(args, "use_original_biophysics", False),
    )
    t, current, soma_v = simulate_mod_current(
        cell_params,
        rows,
        gbar,
        irradiance_scale,
        rs_mohm,
        protocol,
        args.dt_ms,
        1200.0,
        1000.0,
        q10_scales=q10_scales,
        use_photon_flux=use_photon_flux,
        use_original_biophysics=getattr(args, "use_original_biophysics", False),
        subtract_no_light=getattr(args, "use_original_biophysics", False),
    )
    peak, _, time_to_peak, tau = current_metrics(t, current, protocol)

    intensity_values = WANG_INTENSITY_VALUES_MW_MM2
    intensity_peaks = []
    for irr in intensity_values:
        _, current_i, _ = simulate_mod_current(
            cell_params,
            rows,
            gbar,
            irradiance_scale,
            rs_mohm,
            protocol,
            args.dt_ms,
            180.0,
            100.0,
            protocol_irradiance_scale=irr / protocol.irradiance_mw_mm2,
            q10_scales=q10_scales,
            use_photon_flux=use_photon_flux,
            use_original_biophysics=getattr(args, "use_original_biophysics", False),
            subtract_no_light=getattr(args, "use_original_biophysics", False),
        )
        intensity_peaks.append(float(np.max(current_i)))
    intensity_k, intensity_imax, intensity_hill_n = fit_hill_k(intensity_values, intensity_peaks)

    duration_values = [1, 2, 3, 4, 5, 8, 10, 20, 50, 100]
    duration_peaks = []
    for dur in duration_values:
        _, current_d, _ = simulate_mod_current(
            cell_params,
            rows,
            gbar,
            irradiance_scale,
            rs_mohm,
            protocol,
            args.dt_ms,
            max(120.0, float(dur) + 60.0),
            float(dur),
            q10_scales=q10_scales,
            use_photon_flux=use_photon_flux,
            use_original_biophysics=getattr(args, "use_original_biophysics", False),
            subtract_no_light=getattr(args, "use_original_biophysics", False),
        )
        duration_peaks.append(float(np.max(current_d)))

    total_area = sum(row["area_um2"] for row in rows)
    weighted_area = sum(row["area_um2"] * row.get("localization_weight", 1.0) for row in rows)
    mean_irr = sum(row["irradiance_mw_mm2"] * row["area_um2"] for row in rows) / total_area
    result = ModSuiteResult(
        label=label,
        localization=args.localization,
        calibration_target=args.calibration_target,
        calibration_irradiance_mw_mm2=float(calibration_irradiance_mw_mm2),
        calibration_duration_ms=float(calibration_duration_ms),
        calibration_target_current_nA=float(target_current_nA),
        rs_mohm=rs_mohm,
        irradiance_scale=irradiance_scale,
        q10=q10,
        q10_rate_scale=q10_rate_scale,
        q10_mode=args.q10_mode,
        gbar_mS_cm2=float(gbar),
        peak_1s_9p2_nA=float(peak),
        time_to_peak_1s_9p2_ms=float(time_to_peak),
        inactivation_tau_1s_9p2_ms=tau,
        soma_v_min_mV=float(np.min(soma_v)),
        soma_v_max_mV=float(np.max(soma_v)),
        soma_v_mean_mV=float(np.mean(soma_v)),
        intensity_k_mw_mm2=intensity_k,
        intensity_imax_nA=float(intensity_imax),
        intensity_hill_n=float(intensity_hill_n),
        duration_k_ms=half_max_x([float(v) for v in duration_values], duration_peaks),
        duration_imax_nA=float(max(duration_peaks)),
        retained_area_um2=float(total_area),
        weighted_area_um2=float(weighted_area),
        mean_irradiance_mw_mm2=float(mean_irr),
    )
    return result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks


def plot_condition(out, result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks, protocol):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    axes[0].plot(t, current, color="#222222", lw=1.3)
    axes[0].axhline(protocol.target_peak_current_nA, color="#d62728", ls="--", lw=1)
    axes[0].set_xlim(-10, 250)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("SEClamp current magnitude (nA)")
    axes[0].set_title(f"{result.label}, Rs={result.rs_mohm:g}, scale={result.irradiance_scale:g}")

    axes[1].plot(intensity_values, intensity_peaks, marker="o")
    axes[1].axvline(0.84, color="#d62728", ls="--", lw=1)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("irradiance (mW/mm2)")
    axes[1].set_ylabel("peak current (nA)")
    axes[1].set_title(f"intensity K={result.intensity_k_mw_mm2:.3g}")

    axes[2].plot(duration_values, duration_peaks, marker="o")
    axes[2].axvline(3.2, color="#d62728", ls="--", lw=1)
    axes[2].set_xscale("log")
    axes[2].set_xlabel("duration (ms)")
    axes[2].set_ylabel("peak current (nA)")
    axes[2].set_title(f"duration K={result.duration_k_ms:.3g}")
    name = f"{result.label}_{result.localization}_Rs{str(result.rs_mohm).replace('.', 'p')}_scale{str(result.irradiance_scale).replace('.', 'p')}_mod_suite.png"
    fig.savefig(out / name, dpi=220)
    plt.close(fig)


def plot_localization_summary(results, out: Path):
    if not results:
        return
    labels = [r.localization for r in results]
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.8), constrained_layout=True)
    metrics = [
        ("gbar_mS_cm2", "gbar (mS/cm2)"),
        ("time_to_peak_1s_9p2_ms", "time to peak (ms)"),
        ("intensity_k_mw_mm2", "intensity K (mW/mm2)"),
        ("duration_k_ms", "duration K (ms)"),
    ]
    for ax, (field, ylabel) in zip(axes, metrics):
        values = [getattr(r, field) for r in results]
        ax.bar(labels, values, color="#4c78a8")
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=35)
    fig.savefig(out / "localization_summary.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang2007_seclamp_mod_suite"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--pt-cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--it-cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/IT5B_full_cellParams.pkl"))
    parser.add_argument("--cells", choices=["pt", "it", "both"], default="both")
    parser.add_argument("--rs-mohm", type=float, default=5.0)
    parser.add_argument("--irradiance-scale", type=float, default=1.0)
    parser.add_argument("--q10-values", default="1.0")
    parser.add_argument("--q10-mode", choices=["generic", "williams"], default="generic")
    parser.add_argument("--activation-mode", choices=["linear", "photon"], default="linear")
    parser.add_argument("--localization", default="uniform")
    parser.add_argument("--localization-suite", action="store_true")
    parser.add_argument("--proximal-cutoff-um", type=float, default=150.0)
    parser.add_argument("--calibration-target", choices=["peak_9p2_1s", "imax_saturating"], default="peak_9p2_1s")
    parser.add_argument("--saturating-irradiance-mw-mm2", type=float, default=100.0)
    parser.add_argument("--saturating-duration-ms", type=float, default=100.0)
    parser.add_argument("--reference-c", type=float, default=37.0)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--soma-entry-depth-um", type=float, default=100.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--illumination-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--binary-iterations", type=int, default=12)
    parser.add_argument("--use-original-biophysics", action="store_true")
    args = parser.parse_args()

    load_mechanisms(args.repo_root.resolve())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    protocol = Wang2007Protocol()
    cell_items = []
    if args.cells in {"pt", "both"}:
        cell_items.append(("PT5B_full", args.pt_cell_params))
    if args.cells in {"it", "both"}:
        cell_items.append(("IT5B_full", args.it_cell_params))

    results = []
    localization_modes = (
        ["uniform", "soma-only", "dendrite-only", "proximal-enriched", "proximal-reduced"]
        if args.localization_suite
        else [args.localization]
    )
    for label, path in cell_items:
        h.load_file("stdrun.hoc")
        template_cell = CellFromNetPyNE(path)
        rows = prepare_template_rows(template_cell, protocol, args)
        for localization in localization_modes:
            args.localization = localization
            for q10 in [float(raw.strip()) for raw in args.q10_values.split(",") if raw.strip()]:
                result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks = run_condition(
                    label, path, rows, args.rs_mohm, args.irradiance_scale, q10, protocol, args
                )
                results.append(result)
                plot_condition(args.out_dir, result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks, protocol)
                print(result)
    if args.localization_suite and len(cell_items) == 1:
        plot_localization_summary(results, args.out_dir)

    with (args.out_dir / "seclamp_mod_suite_results.json").open("w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)
    with (args.out_dir / "seclamp_mod_suite_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in results)


if __name__ == "__main__":
    main()
