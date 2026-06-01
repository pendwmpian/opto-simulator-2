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

from run_pt5b_current_injection import simulate_current_step
from run_wang_seclamp_mod_suite import load_mechanisms


def parse_floats(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def summarize_trace(trace, test_current_nA: float, measure_start_ms: float, measure_stop_ms: float, step_duration_ms: float, threshold_mV: float):
    t = trace["t"]
    v = trace["v"]
    baseline = v[(t >= 0.0) & (t <= min(20.0, measure_start_ms))]
    measured = v[(t >= measure_start_ms) & (t <= measure_stop_ms)]
    rest_v = float(np.mean(baseline)) if len(baseline) else float(v[0])
    steady_v = float(np.mean(measured)) if len(measured) else float("nan")
    above = v >= threshold_mV
    spike_count = 0
    was_above = False
    for time, is_above in zip(t, above):
        if time > step_duration_ms:
            break
        if is_above and not was_above:
            spike_count += 1
        was_above = bool(is_above)
    return rest_v, steady_v, (steady_v - rest_v) / test_current_nA if test_current_nA else float("nan"), spike_count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/leak_conductance_input_resistance"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--g-pas-scales", default="1.0,0.8,0.6,0.5,0.4,0.3,0.2")
    parser.add_argument("--fi-currents-nA", default="0,0.1,0.2,0.3,0.4,0.5")
    parser.add_argument("--test-current-nA", type=float, default=0.05)
    parser.add_argument("--soma-entry-depth-um", type=float, default=25.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--initial-v-mV", type=float, default=-70.0)
    parser.add_argument("--equilibration-ms", type=float, default=500.0)
    parser.add_argument("--step-duration-ms", type=float, default=300.0)
    parser.add_argument("--measure-start-ms", type=float, default=200.0)
    parser.add_argument("--measure-stop-ms", type=float, default=300.0)
    parser.add_argument("--fi-step-duration-ms", type=float, default=1000.0)
    parser.add_argument("--post-ms", type=float, default=50.0)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())

    scales = parse_floats(args.g_pas_scales)
    fi_currents = parse_floats(args.fi_currents_nA)
    rows = []
    fi_rows = []
    for scale in scales:
        rin_args = argparse.Namespace(**vars(args))
        rin_args.use_original_biophysics = True
        rin_args.apply_slice_pruning = True
        rin_args.g_pas_scale = scale
        rin_args.e_pas_mV = None
        rin_args.holding_current_nA = 0.0
        rin_args.refractory_ms = 2.0
        rin_args.current_values_nA = str(args.test_current_nA)
        trace0 = simulate_current_step(args.cell_params, 0.0, rin_args)
        trace_i = simulate_current_step(args.cell_params, args.test_current_nA, rin_args)
        rest0 = float(np.mean(trace0["v"][(trace0["t"] >= 0.0) & (trace0["t"] <= 20.0)]))
        _, steady_i, _, spike_count = summarize_trace(
            trace_i,
            args.test_current_nA,
            args.measure_start_ms,
            args.measure_stop_ms,
            args.step_duration_ms,
            args.spike_threshold_mV,
        )
        rin = (steady_i - rest0) / args.test_current_nA
        rows.append(
            {
                "label": args.label,
                "g_pas_scale": scale,
                "soma_entry_depth_um": args.soma_entry_depth_um,
                "test_current_nA": args.test_current_nA,
                "rest_v_mV": rest0,
                "steady_v_mV": steady_i,
                "delta_v_mV": steady_i - rest0,
                "input_resistance_MOhm": rin,
                "spike_count_test_step": spike_count,
                "retained_area_um2": trace_i["retained_area_um2"],
                "retained_area_fraction": trace_i["retained_area_fraction"],
            }
        )
        fi_args = argparse.Namespace(**vars(rin_args))
        fi_args.step_duration_ms = args.fi_step_duration_ms
        for current in fi_currents:
            trace = simulate_current_step(args.cell_params, current, fi_args)
            rest_v, _, _, spikes = summarize_trace(
                trace,
                current if current else 1.0,
                min(200.0, args.fi_step_duration_ms / 2.0),
                args.fi_step_duration_ms,
                args.fi_step_duration_ms,
                args.spike_threshold_mV,
            )
            fi_rows.append(
                {
                    "label": args.label,
                    "g_pas_scale": scale,
                    "current_nA": current,
                    "rest_v_mV": rest_v,
                    "spike_count_1s": spikes,
                    "firing_rate_Hz": spikes / (args.fi_step_duration_ms / 1000.0),
                }
            )

    with (args.out_dir / "leak_conductance_input_resistance.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with (args.out_dir / "leak_conductance_fi.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fi_rows[0].keys()))
        writer.writeheader()
        writer.writerows(fi_rows)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), constrained_layout=True)
    axes[0].plot([row["g_pas_scale"] for row in rows], [row["input_resistance_MOhm"] for row in rows], marker="o")
    axes[0].axhline(126.0, color="#d62728", ls="--", lw=1, label="Wang ChR2+")
    axes[0].axhline(110.0, color="#999999", ls=":", lw=1, label="Wang ChR2-")
    axes[0].invert_xaxis()
    axes[0].set_xlabel("g_pas scale")
    axes[0].set_ylabel("Rin (MOhm)")
    axes[0].set_title("Input resistance")
    axes[0].legend(fontsize=8)

    axes[1].plot([row["g_pas_scale"] for row in rows], [row["rest_v_mV"] for row in rows], marker="o")
    axes[1].axhline(-58.3, color="#d62728", ls="--", lw=1)
    axes[1].invert_xaxis()
    axes[1].set_xlabel("g_pas scale")
    axes[1].set_ylabel("rest Vm (mV)")
    axes[1].set_title("Resting voltage")

    for current in fi_currents:
        sub = [row for row in fi_rows if row["current_nA"] == current]
        axes[2].plot([row["g_pas_scale"] for row in sub], [row["spike_count_1s"] for row in sub], marker="o", label=f"{current:g} nA")
    axes[2].invert_xaxis()
    axes[2].set_xlabel("g_pas scale")
    axes[2].set_ylabel("spikes / 1 s")
    axes[2].set_title("f-I check")
    axes[2].legend(fontsize=7, ncol=2)
    fig.savefig(args.out_dir / "leak_conductance_input_resistance.png", dpi=220)
    plt.close(fig)

    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
