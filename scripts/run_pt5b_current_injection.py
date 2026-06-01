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


def prune_to_slice(cell: CellFromNetPyNE, args) -> tuple[float, float]:
    protocol = Wang2007Protocol()
    soma_axis = soma_axis_value(cell, args.slice_axis)
    slice_center_um = soma_axis + protocol.slice_thickness_um / 2.0 - args.soma_entry_depth_um
    rows, _ = segment_rows_vertical_slice(
        cell,
        args.slice_axis,
        slice_center_um,
        protocol.slice_thickness_um,
        args.truncate_radius_um,
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


def adjust_passive(cell: CellFromNetPyNE, g_pas_scale: float, e_pas_mV: float | None) -> None:
    for sec in cell.sections.values():
        if not h.ismembrane("pas", sec=sec):
            continue
        sec.g_pas *= g_pas_scale
        if e_pas_mV is not None:
            sec.e_pas = e_pas_mV


def spike_peaks(t: np.ndarray, v: np.ndarray, threshold_mV: float, refractory_ms: float):
    peaks = []
    last_t = -1.0e9
    for i in range(1, len(v) - 1):
        if v[i] < threshold_mV:
            continue
        if v[i] >= v[i - 1] and v[i] > v[i + 1] and t[i] - last_t >= refractory_ms:
            peaks.append((float(t[i]), float(v[i])))
            last_t = float(t[i])
    return peaks


def simulate_current_step(cell_params: Path, amp_nA: float, args):
    h.load_file("stdrun.hoc")
    h.celsius = args.target_c
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=args.use_original_biophysics)
    retained_area_um2 = float("nan")
    retained_area_fraction = float("nan")
    if args.apply_slice_pruning:
        retained_area_um2, retained_area_fraction = prune_to_slice(cell, args)
    adjust_passive(cell, args.g_pas_scale, args.e_pas_mV)
    soma = soma_section(cell)(0.5)

    hold = h.IClamp(soma)
    hold.delay = -1.0e9
    hold.dur = 1.0e9
    hold.amp = args.holding_current_nA

    stim = h.IClamp(soma)
    stim.delay = 0.0
    stim.dur = args.step_duration_ms
    stim.amp = amp_nA

    h.dt = args.dt_ms
    h.tstop = args.post_ms
    h.finitialize(args.initial_v_mV)
    h.t = -args.equilibration_ms

    times = []
    voltages = []
    currents = []
    while h.t <= args.step_duration_ms + args.post_ms:
        if h.t >= 0.0:
            times.append(float(h.t))
            voltages.append(float(soma.v))
            currents.append(args.holding_current_nA + (amp_nA if h.t < args.step_duration_ms else 0.0))
        h.fadvance()

    return {
        "t": np.asarray(times),
        "v": np.asarray(voltages),
        "i": np.asarray(currents),
        "retained_area_um2": retained_area_um2,
        "retained_area_fraction": retained_area_fraction,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/pt5b_current_injection"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--current-values-nA", default="0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5,2.0")
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--initial-v-mV", type=float, default=-70.0)
    parser.add_argument("--equilibration-ms", type=float, default=500.0)
    parser.add_argument("--step-duration-ms", type=float, default=1000.0)
    parser.add_argument("--post-ms", type=float, default=200.0)
    parser.add_argument("--dt-ms", type=float, default=0.05)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    parser.add_argument("--refractory-ms", type=float, default=2.0)
    parser.add_argument("--use-original-biophysics", action="store_true")
    parser.add_argument("--g-pas-scale", type=float, default=1.0)
    parser.add_argument("--e-pas-mV", type=float, default=None)
    parser.add_argument("--holding-current-nA", type=float, default=0.0)
    parser.add_argument("--apply-slice-pruning", action="store_true")
    parser.add_argument("--soma-entry-depth-um", type=float, default=100.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    currents = [float(raw.strip()) for raw in args.current_values_nA.split(",") if raw.strip()]

    rows = []
    traces = {}
    for amp in currents:
        trace = simulate_current_step(args.cell_params, amp, args)
        spikes = spike_peaks(trace["t"], trace["v"], args.spike_threshold_mV, args.refractory_ms)
        during = trace["t"] <= args.step_duration_ms
        spike_times = [spike_t for spike_t, _ in spikes if spike_t <= args.step_duration_ms]
        rows.append(
            {
                "label": args.label,
                "current_nA": amp,
                "rest_v_mV": float(trace["v"][0]),
                "max_v_mV": float(np.max(trace["v"][during])),
                "min_v_mV": float(np.min(trace["v"][during])),
                "spike_count_1s": len(spike_times),
                "firing_rate_Hz": len(spike_times) / (args.step_duration_ms / 1000.0),
                "first_spike_latency_ms": spike_times[0] if spike_times else float("nan"),
                "last_spike_time_ms": spike_times[-1] if spike_times else float("nan"),
                "retained_area_um2": trace["retained_area_um2"],
                "retained_area_fraction": trace["retained_area_fraction"],
            }
        )
        traces[amp] = trace

    baseline_rows = [row for row in rows if abs(row["current_nA"]) < 1.0e-12]
    baseline_max_v = baseline_rows[0]["max_v_mV"] if baseline_rows else rows[0]["max_v_mV"]
    for row in rows:
        if abs(row["current_nA"]) < 1.0e-12:
            row["apparent_input_resistance_MOhm"] = float("nan")
        else:
            row["apparent_input_resistance_MOhm"] = (row["max_v_mV"] - baseline_max_v) / row["current_nA"]

    with (args.out_dir / "current_injection_fi.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), constrained_layout=True)
    example_currents = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
    for amp in example_currents:
        if amp in traces:
            trace = traces[amp]
            axes[0].plot(trace["t"], trace["v"], lw=1.0, label=f"{amp:g} nA")
    axes[0].set_xlim(-10, args.step_duration_ms + 50.0)
    axes[0].set_xlabel("time (ms)")
    axes[0].set_ylabel("soma Vm (mV)")
    axes[0].set_title("Somatic current steps")
    axes[0].legend(fontsize=8)

    axes[1].plot([row["current_nA"] for row in rows], [row["spike_count_1s"] for row in rows], marker="o")
    axes[1].set_xlabel("current (nA)")
    axes[1].set_ylabel("spikes / 1 s")
    axes[1].set_title("f-I count")

    axes[2].plot([row["current_nA"] for row in rows], [row["first_spike_latency_ms"] for row in rows], marker="o")
    axes[2].set_xlabel("current (nA)")
    axes[2].set_ylabel("first spike latency (ms)")
    axes[2].set_title("latency")
    fig.savefig(args.out_dir / "current_injection_fi.png", dpi=220)
    plt.close(fig)

    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
