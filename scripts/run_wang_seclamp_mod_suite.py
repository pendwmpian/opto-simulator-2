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


@dataclass(frozen=True)
class ModSuiteResult:
    label: str
    rs_mohm: float
    irradiance_scale: float
    gbar_mS_cm2: float
    peak_1s_9p2_nA: float
    time_to_peak_1s_9p2_ms: float
    inactivation_tau_1s_9p2_ms: float | None
    intensity_k_mw_mm2: float | None
    intensity_imax_nA: float
    duration_k_ms: float | None
    duration_imax_nA: float
    retained_area_um2: float
    mean_irradiance_mw_mm2: float


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


def install_chr2(cell: CellFromNetPyNE, rows_template, gbar_mS_cm2: float, irradiance_scale: float):
    for sec in cell.sections.values():
        sec.insert("chr2_4state")
        for seg in sec:
            seg.chr2_4state.gbar = 0.0
            seg.chr2_4state.irr = 0.0

    rows = []
    for template in rows_template:
        sec_name = template["sec"].name().split(".")[-1]
        sec = cell.sections[sec_name]
        seg = sec(template["xloc"])
        seg.chr2_4state.gbar = gbar_mS_cm2
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
):
    h.load_file("stdrun.hoc")
    cell = CellFromNetPyNE(cell_params)
    rows = install_chr2(cell, rows_template, gbar_mS_cm2, irradiance_scale)

    clamp = h.SEClamp(soma_section(cell)(0.5))
    clamp.dur1 = tstop_ms + 1.0
    clamp.amp1 = protocol.holding_voltage_mV
    clamp.rs = rs_mohm

    h.dt = dt_ms
    h.tstop = tstop_ms
    h.finitialize(protocol.holding_voltage_mV)
    times = []
    currents = []
    while h.t <= tstop_ms:
        light_on = 0.0 <= h.t < pulse_duration_ms
        for row in rows:
            row["seg"].chr2_4state.irr = row["base_irradiance_mw_mm2"] * protocol_irradiance_scale if light_on else 0.0
        times.append(float(h.t))
        currents.append(float(abs(clamp.i)))
        h.fadvance()
    return np.asarray(times), np.asarray(currents)


def calibrate_gbar(cell_params, rows, irradiance_scale, rs_mohm, protocol, dt_ms, iterations):
    lo = 1.0e-5
    hi = 1.0
    last = None
    for _ in range(iterations):
        mid = (lo * hi) ** 0.5
        t, current = simulate_mod_current(
            cell_params,
            rows,
            mid,
            irradiance_scale,
            rs_mohm,
            protocol,
            dt_ms,
            1200.0,
            1000.0,
        )
        peak = float(np.max(current))
        last = (mid, t, current)
        if peak < protocol.target_peak_current_nA:
            lo = mid
        else:
            hi = mid
    if last is None:
        raise RuntimeError("Calibration did not run")
    return last


def run_condition(label, cell_params, rows, rs_mohm, irradiance_scale, protocol, args):
    gbar, t, current = calibrate_gbar(cell_params, rows, irradiance_scale, rs_mohm, protocol, args.dt_ms, args.binary_iterations)
    peak, _, time_to_peak, tau = current_metrics(t, current, protocol)

    intensity_values = [0.05, 0.1, 0.2, 0.4, 0.84, 1.0, 2.0, 4.6, 9.2]
    intensity_peaks = []
    for irr in intensity_values:
        _, current_i = simulate_mod_current(
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
        )
        intensity_peaks.append(float(np.max(current_i)))

    duration_values = [1, 2, 3, 4, 5, 8, 10, 20, 50, 100]
    duration_peaks = []
    for dur in duration_values:
        _, current_d = simulate_mod_current(
            cell_params,
            rows,
            gbar,
            irradiance_scale,
            rs_mohm,
            protocol,
            args.dt_ms,
            max(120.0, float(dur) + 60.0),
            float(dur),
        )
        duration_peaks.append(float(np.max(current_d)))

    total_area = sum(row["area_um2"] for row in rows)
    mean_irr = sum(row["irradiance_mw_mm2"] * row["area_um2"] for row in rows) / total_area
    result = ModSuiteResult(
        label=label,
        rs_mohm=rs_mohm,
        irradiance_scale=irradiance_scale,
        gbar_mS_cm2=float(gbar),
        peak_1s_9p2_nA=float(peak),
        time_to_peak_1s_9p2_ms=float(time_to_peak),
        inactivation_tau_1s_9p2_ms=tau,
        intensity_k_mw_mm2=half_max_x(intensity_values, intensity_peaks),
        intensity_imax_nA=float(max(intensity_peaks)),
        duration_k_ms=half_max_x([float(v) for v in duration_values], duration_peaks),
        duration_imax_nA=float(max(duration_peaks)),
        retained_area_um2=float(total_area),
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
    name = f"{result.label}_Rs{str(result.rs_mohm).replace('.', 'p')}_scale{str(result.irradiance_scale).replace('.', 'p')}_mod_suite.png"
    fig.savefig(out / name, dpi=220)
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
    parser.add_argument("--soma-entry-depth-um", type=float, default=100.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--illumination-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--vitro-mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--binary-iterations", type=int, default=12)
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
    for label, path in cell_items:
        h.load_file("stdrun.hoc")
        template_cell = CellFromNetPyNE(path)
        rows = prepare_template_rows(template_cell, protocol, args)
        result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks = run_condition(
            label, path, rows, args.rs_mohm, args.irradiance_scale, protocol, args
        )
        results.append(result)
        plot_condition(args.out_dir, result, t, current, intensity_values, intensity_peaks, duration_values, duration_peaks, protocol)
        print(result)

    with (args.out_dir / "seclamp_mod_suite_results.json").open("w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)
    with (args.out_dir / "seclamp_mod_suite_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in results)


if __name__ == "__main__":
    main()
