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

from calibrate_wang2007_gbar import Wang2007Protocol, segment_rows_vertical_slice, soma_axis_value
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_seclamp_mod_suite import load_mechanisms, soma_section


def electrically_prune_to_slice(cell: CellFromNetPyNE, slice_axis: str, soma_entry_depth_um: float, truncate_radius_um: float | None):
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
    kept = {(row["sec"].name().split(".")[-1], round(float(row["xloc"]), 12)) for row in rows}
    kept_area = 0.0
    total_area = 0.0
    for sec_name, sec in cell.sections.items():
        for seg in sec:
            area_um2 = float(h.area(seg.x, sec=sec))
            total_area += area_um2
            key = (sec_name, round(float(seg.x), 12))
            if key in kept:
                kept_area += area_um2
                continue
            for mech_name in list(sec.psection().get("density_mechs", {}).keys()):
                mech = getattr(seg, mech_name)
                for attr in ("g", "gbar", "gcalbar", "gcanbar", "gmax"):
                    if hasattr(mech, attr):
                        try:
                            setattr(mech, attr, 0.0)
                        except Exception:
                            pass
            sec.cm = 1.0e-9
    return kept_area, kept_area / total_area if total_area else 0.0


def measure_rin(cell_params: Path, args, soma_entry_depth_um: float | None):
    h.load_file("stdrun.hoc")
    h.celsius = args.target_c
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=True)

    retained_area_um2 = float("nan")
    retained_area_fraction = float("nan")
    if soma_entry_depth_um is not None:
        retained_area_um2, retained_area_fraction = electrically_prune_to_slice(
            cell,
            args.slice_axis,
            soma_entry_depth_um,
            args.truncate_radius_um,
        )

    soma = soma_section(cell)(0.5)
    stim = h.IClamp(soma)
    stim.delay = 0.0
    stim.dur = args.step_duration_ms
    stim.amp = args.test_current_nA

    h.dt = args.dt_ms
    h.tstop = args.step_duration_ms + args.post_ms
    h.finitialize(args.initial_v_mV)
    h.t = -args.equilibration_ms

    baseline_v = []
    step_v = []
    spike_count = 0
    above = False
    while h.t <= h.tstop:
        if -100.0 <= h.t < 0.0:
            baseline_v.append(float(soma.v))
        if args.measure_start_ms <= h.t <= args.measure_stop_ms:
            step_v.append(float(soma.v))
        if 0.0 <= h.t <= args.step_duration_ms:
            now_above = float(soma.v) >= args.spike_threshold_mV
            if now_above and not above:
                spike_count += 1
            above = now_above
        h.fadvance()

    rest_v = float(np.mean(baseline_v)) if baseline_v else float(soma.v)
    steady_v = float(np.mean(step_v)) if step_v else float("nan")
    rin = (steady_v - rest_v) / args.test_current_nA if args.test_current_nA else float("nan")
    return {
        "soma_entry_depth_um": "full" if soma_entry_depth_um is None else soma_entry_depth_um,
        "slice_axis": args.slice_axis if soma_entry_depth_um is not None else "none",
        "test_current_nA": args.test_current_nA,
        "rest_v_mV": rest_v,
        "steady_v_mV": steady_v,
        "delta_v_mV": steady_v - rest_v,
        "input_resistance_MOhm": rin,
        "spike_count": spike_count,
        "retained_area_um2": retained_area_um2,
        "retained_area_fraction": retained_area_fraction,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/slice_position_input_resistance"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--soma-entry-depths-um", default="25,50,75,100,125,150,200,250")
    parser.add_argument("--include-full", action="store_true")
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--initial-v-mV", type=float, default=-70.0)
    parser.add_argument("--equilibration-ms", type=float, default=500.0)
    parser.add_argument("--test-current-nA", type=float, default=0.05)
    parser.add_argument("--step-duration-ms", type=float, default=300.0)
    parser.add_argument("--measure-start-ms", type=float, default=200.0)
    parser.add_argument("--measure-stop-ms", type=float, default=300.0)
    parser.add_argument("--post-ms", type=float, default=50.0)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())

    depths: list[float | None] = [float(raw.strip()) for raw in args.soma_entry_depths_um.split(",") if raw.strip()]
    if args.include_full:
        depths = [None] + depths

    rows = [measure_rin(args.cell_params, args, depth) for depth in depths]

    csv_path = args.out_dir / "slice_position_input_resistance.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    plot_rows = [row for row in rows if row["soma_entry_depth_um"] != "full"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), constrained_layout=True)
    xs = [float(row["soma_entry_depth_um"]) for row in plot_rows]
    axes[0].plot(xs, [row["input_resistance_MOhm"] for row in plot_rows], marker="o")
    axes[0].axhline(126.0, color="#d62728", ls="--", lw=1, label="Wang ChR2+")
    axes[0].axhline(110.0, color="#999999", ls=":", lw=1, label="Wang ChR2-")
    axes[0].set_xlabel("soma entry depth (um)")
    axes[0].set_ylabel("Rin (MOhm)")
    axes[0].set_title("Input resistance")
    axes[0].legend(fontsize=8)

    axes[1].plot(xs, [row["retained_area_um2"] for row in plot_rows], marker="o")
    axes[1].set_xlabel("soma entry depth (um)")
    axes[1].set_ylabel("retained area (um2)")
    axes[1].set_title("Retained membrane area")

    axes[2].plot(xs, [row["rest_v_mV"] for row in plot_rows], marker="o")
    axes[2].axhline(-58.3, color="#d62728", ls="--", lw=1)
    axes[2].set_xlabel("soma entry depth (um)")
    axes[2].set_ylabel("rest Vm (mV)")
    axes[2].set_title("Resting voltage")
    fig.savefig(args.out_dir / "slice_position_input_resistance.png", dpi=220)
    plt.close(fig)

    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
