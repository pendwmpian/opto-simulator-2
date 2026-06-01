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

from calibrate_wang2007_gbar import Wang2007Protocol, segment_rows_vertical_slice, soma_axis_value
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_seclamp_mod_suite import load_mechanisms, soma_section


@dataclass(frozen=True)
class RinResult:
    label: str
    slice_axis: str
    soma_entry_depth_um: float
    slice_center_um: float
    truncate_radius_um: float
    current_nA: float
    baseline_v_mV: float
    steady_v_mV: float
    delta_v_mV: float
    input_resistance_MOhm: float
    target_input_resistance_MOhm: float
    abs_error_MOhm: float
    retained_area_um2: float
    total_area_um2: float
    retained_area_fraction: float


def kept_segment_keys(cell: CellFromNetPyNE, slice_axis: str, soma_entry_depth_um: float, truncate_radius_um: float):
    protocol = Wang2007Protocol()
    soma_axis = soma_axis_value(cell, slice_axis)
    slice_center_um = soma_axis + protocol.slice_thickness_um / 2.0 - soma_entry_depth_um
    rows, _ = segment_rows_vertical_slice(
        cell,
        slice_axis,
        slice_center_um,
        protocol.slice_thickness_um,
        truncate_radius_um,
    )
    keys = {(row["sec"].name().split(".")[-1], round(float(row["xloc"]), 12)) for row in rows}
    return keys, slice_center_um


def suppress_offslice_membrane(cell: CellFromNetPyNE, kept_keys: set[tuple[str, float]]) -> tuple[float, float]:
    retained_area_um2 = 0.0
    total_area_um2 = 0.0
    for sec_name, sec in cell.sections.items():
        sec_retained = False
        for seg in sec:
            area_um2 = float(h.area(seg.x, sec=sec))
            total_area_um2 += area_um2
            key = (sec_name, round(float(seg.x), 12))
            if key in kept_keys:
                retained_area_um2 += area_um2
                sec_retained = True
                continue
            for mech_name in list(sec.psection().get("density_mechs", {}).keys()):
                mech = getattr(seg, mech_name)
                for attr in ("g", "gbar", "gcalbar", "gcanbar", "gmax"):
                    if hasattr(mech, attr):
                        try:
                            setattr(mech, attr, 0.0)
                        except Exception:
                            pass
        if not sec_retained:
            sec.cm = 1.0e-9
    return retained_area_um2, total_area_um2


def scale_passive(cell: CellFromNetPyNE, g_pas_scale: float, e_pas_mV: float | None) -> None:
    for sec in cell.sections.values():
        if h.ismembrane("pas", sec=sec):
            sec.g_pas *= g_pas_scale
            if e_pas_mV is not None:
                sec.e_pas = e_pas_mV
        if h.ismembrane("hh", sec=sec):
            sec.gl_hh *= g_pas_scale
            if e_pas_mV is not None:
                sec.el_hh = e_pas_mV


def measure_rin(cell_params: Path, slice_axis: str, soma_entry_depth_um: float, args) -> tuple[RinResult, dict[str, np.ndarray]]:
    h.load_file("stdrun.hoc")
    h.celsius = args.target_c
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=args.use_original_biophysics)
    kept_keys, slice_center_um = kept_segment_keys(cell, slice_axis, soma_entry_depth_um, args.truncate_radius_um)
    retained_area_um2, total_area_um2 = suppress_offslice_membrane(cell, kept_keys)
    scale_passive(cell, args.g_pas_scale, args.e_pas_mV)

    soma = soma_section(cell)(0.5)
    stim = h.IClamp(soma)
    stim.delay = args.step_start_ms
    stim.dur = args.step_duration_ms
    stim.amp = args.current_nA

    h.dt = args.dt_ms
    h.tstop = args.step_start_ms + args.step_duration_ms + args.post_ms
    h.finitialize(args.initial_v_mV)
    h.t = -args.equilibration_ms

    times = []
    voltages = []
    while h.t <= h.tstop:
        if h.t >= 0.0:
            times.append(float(h.t))
            voltages.append(float(soma.v))
        h.fadvance()

    t = np.asarray(times)
    v = np.asarray(voltages)
    baseline_mask = (t >= args.baseline_start_ms) & (t < args.baseline_stop_ms)
    steady_start = args.step_start_ms + args.step_duration_ms - args.steady_window_ms
    steady_mask = (t >= steady_start) & (t < args.step_start_ms + args.step_duration_ms)
    baseline_v = float(np.mean(v[baseline_mask]))
    steady_v = float(np.mean(v[steady_mask]))
    delta_v = steady_v - baseline_v
    rin = delta_v / args.current_nA
    result = RinResult(
        label=args.label,
        slice_axis=slice_axis,
        soma_entry_depth_um=float(soma_entry_depth_um),
        slice_center_um=float(slice_center_um),
        truncate_radius_um=float(args.truncate_radius_um),
        current_nA=float(args.current_nA),
        baseline_v_mV=baseline_v,
        steady_v_mV=steady_v,
        delta_v_mV=delta_v,
        input_resistance_MOhm=float(rin),
        target_input_resistance_MOhm=float(args.target_rin_MOhm),
        abs_error_MOhm=abs(float(rin) - args.target_rin_MOhm),
        retained_area_um2=float(retained_area_um2),
        total_area_um2=float(total_area_um2),
        retained_area_fraction=float(retained_area_um2 / total_area_um2) if total_area_um2 else float("nan"),
    )
    return result, {"t": t, "v": v}


def write_csv(path: Path, rows: list[RinResult]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def plot_summary(out_dir: Path, rows: list[RinResult], traces: dict[tuple[str, float], dict[str, np.ndarray]], best: RinResult) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    axes[0].axhline(best.target_input_resistance_MOhm, color="#d62728", ls="--", lw=1, label="126 MOhm target")
    for axis in sorted({row.slice_axis for row in rows}):
        axis_rows = sorted([row for row in rows if row.slice_axis == axis], key=lambda row: row.soma_entry_depth_um)
        axes[0].plot(
            [row.soma_entry_depth_um for row in axis_rows],
            [row.input_resistance_MOhm for row in axis_rows],
            marker="o",
            label=f"{axis}-thickness",
        )
    axes[0].set_xlabel("soma entry depth (um)")
    axes[0].set_ylabel("input resistance (MOhm)")
    axes[0].set_title("Slice-position Rin sweep")
    axes[0].legend(frameon=False, fontsize=8)

    for axis in sorted({row.slice_axis for row in rows}):
        axis_rows = sorted([row for row in rows if row.slice_axis == axis], key=lambda row: row.soma_entry_depth_um)
        axes[1].plot(
            [row.soma_entry_depth_um for row in axis_rows],
            [row.retained_area_fraction for row in axis_rows],
            marker="o",
            label=f"{axis}-thickness",
        )
    axes[1].set_xlabel("soma entry depth (um)")
    axes[1].set_ylabel("retained membrane area fraction")
    axes[1].set_title("Retained area")
    axes[1].legend(frameon=False, fontsize=8)

    key = (best.slice_axis, best.soma_entry_depth_um)
    trace = traces[key]
    axes[2].plot(trace["t"], trace["v"], color="#222222", lw=1.2)
    axes[2].axvspan(200.0, 1100.0, color="#e5e5e5", zorder=-1)
    axes[2].set_xlabel("time (ms)")
    axes[2].set_ylabel("soma Vm (mV)")
    axes[2].set_title(f"Best trace: {best.slice_axis}, {best.soma_entry_depth_um:g} um")
    fig.savefig(out_dir / "pt5b_slice_rin_sweep.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/pt5b_slice_rin_sweep"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full_ChR2_negative")
    parser.add_argument("--slice-axes", default="z,x")
    parser.add_argument("--soma-entry-depths-um", default="25,50,75,100,125,150,175,200,225,250,275")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--target-rin-MOhm", type=float, default=126.0)
    parser.add_argument("--current-nA", type=float, default=-0.05)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--initial-v-mV", type=float, default=-70.0)
    parser.add_argument("--equilibration-ms", type=float, default=500.0)
    parser.add_argument("--step-start-ms", type=float, default=200.0)
    parser.add_argument("--step-duration-ms", type=float, default=900.0)
    parser.add_argument("--post-ms", type=float, default=200.0)
    parser.add_argument("--baseline-start-ms", type=float, default=100.0)
    parser.add_argument("--baseline-stop-ms", type=float, default=190.0)
    parser.add_argument("--steady-window-ms", type=float, default=100.0)
    parser.add_argument("--dt-ms", type=float, default=0.05)
    parser.add_argument("--use-original-biophysics", action="store_true")
    parser.add_argument("--g-pas-scale", type=float, default=1.0)
    parser.add_argument("--e-pas-mV", type=float, default=None)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    axes = [raw.strip() for raw in args.slice_axes.split(",") if raw.strip()]
    depths = [float(raw.strip()) for raw in args.soma_entry_depths_um.split(",") if raw.strip()]

    rows = []
    traces = {}
    for axis in axes:
        for depth in depths:
            result, trace = measure_rin(args.cell_params, axis, depth, args)
            rows.append(result)
            traces[(axis, depth)] = trace
            print(asdict(result))

    rows.sort(key=lambda row: (row.abs_error_MOhm, row.slice_axis, row.soma_entry_depth_um))
    best = rows[0]
    write_csv(args.out_dir / "pt5b_slice_rin_sweep.csv", rows)
    with (args.out_dir / "best_pt5b_slice_rin.json").open("w") as f:
        json.dump(asdict(best), f, indent=2)
    plot_summary(args.out_dir, rows, traces, best)
    np.savetxt(
        args.out_dir / "best_trace.csv",
        np.column_stack([traces[(best.slice_axis, best.soma_entry_depth_um)]["t"], traces[(best.slice_axis, best.soma_entry_depth_um)]["v"]]),
        delimiter=",",
        header="time_ms,soma_v_mV",
        comments="",
    )


if __name__ == "__main__":
    main()
