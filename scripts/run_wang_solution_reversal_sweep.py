from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from neuron import h

from run_pt5b_current_injection import prune_to_slice, spike_peaks
from run_chr2_single_cell_sweep import CellFromNetPyNE
from run_wang_seclamp_mod_suite import load_mechanisms, soma_section


@dataclass(frozen=True)
class IonProfile:
    name: str
    ena_mV: float | None
    ek_mV: float | None
    eca_mV: float | None
    note: str


def nernst_mV(out_mM: float, in_mM: float, z: float, temp_c: float) -> float:
    r = 8.31446261815324
    f = 96485.33212
    temp_k = temp_c + 273.15
    return 1000.0 * r * temp_k / (z * f) * math.log(out_mM / in_mM)


def wang_ion_profiles(temp_c: float) -> list[IonProfile]:
    na_out = 125.0 + 26.0 + 1.25
    na_in_total = 2.0 + 2.0 * 4.0 + 2.0 * 0.4
    na_in_free2 = 2.0
    k_out = 2.5
    k_in = 130.0 + 0.25
    ca_out = 2.0
    ca_in_free = 0.0001

    ena_total = nernst_mV(na_out, na_in_total, 1.0, temp_c)
    ena_free2 = nernst_mV(na_out, na_in_free2, 1.0, temp_c)
    ek = nernst_mV(k_out, k_in, 1.0, temp_c)
    eca = nernst_mV(ca_out, ca_in_free, 2.0, temp_c)
    return [
        IonProfile("dura_default", None, None, None, "Original Dura-Bernal ion settings."),
        IonProfile(
            "wang_total_na",
            ena_total,
            ek,
            eca,
            "Na_in includes NaCl + Na2ATP + NaGTP sodium; Ca_in free set to 100 nM.",
        ),
        IonProfile(
            "wang_free_na2",
            ena_free2,
            ek,
            eca,
            "Na_in uses only the 2 mM NaCl term; upper-bound E_Na sensitivity check.",
        ),
        IonProfile(
            "wang_total_na_dura_ca",
            ena_total,
            ek,
            None,
            "Na/K from Wang solution; Ca reversal left at Dura-Bernal effective value.",
        ),
    ]


def apply_ion_profile(cell: CellFromNetPyNE, profile: IonProfile) -> None:
    for sec in cell.sections.values():
        ions = sec.psection().get("ions", {})
        if profile.ena_mV is not None and "na" in ions:
            sec.ena = profile.ena_mV
        if profile.ek_mV is not None and "k" in ions:
            sec.ek = profile.ek_mV
        if profile.eca_mV is not None and "ca" in ions:
            sec.eca = profile.eca_mV


def parse_floats(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def simulate_step(cell_params: Path, current_nA: float, profile: IonProfile, args):
    h.load_file("stdrun.hoc")
    h.celsius = args.target_c
    cell = CellFromNetPyNE(cell_params, use_original_biophysics=True)
    retained_area_um2 = float("nan")
    retained_area_fraction = float("nan")
    if args.apply_slice_pruning:
        retained_area_um2, retained_area_fraction = prune_to_slice(cell, args)
    apply_ion_profile(cell, profile)

    soma = soma_section(cell)(0.5)
    stim = h.IClamp(soma)
    stim.delay = 0.0
    stim.dur = args.step_duration_ms
    stim.amp = current_nA

    h.dt = args.dt_ms
    h.tstop = args.step_duration_ms + args.post_ms
    h.finitialize(args.initial_v_mV)
    h.t = -args.equilibration_ms

    times = []
    voltages = []
    while h.t <= h.tstop:
        if h.t >= 0.0:
            times.append(float(h.t))
            voltages.append(float(soma.v))
        h.fadvance()

    return {
        "t": np.asarray(times),
        "v": np.asarray(voltages),
        "retained_area_um2": retained_area_um2,
        "retained_area_fraction": retained_area_fraction,
    }


def baseline_mean(trace, start_ms: float, stop_ms: float) -> float:
    mask = (trace["t"] >= start_ms) & (trace["t"] <= stop_ms)
    if np.any(mask):
        return float(np.mean(trace["v"][mask]))
    return float(trace["v"][0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/wang_solution_reversal_sweep"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--label", default="PT5B_full")
    parser.add_argument("--current-values-nA", default="0,0.1,0.2,0.3,0.4,0.5,0.75")
    parser.add_argument("--test-current-nA", type=float, default=0.05)
    parser.add_argument("--target-c", type=float, default=22.0)
    parser.add_argument("--initial-v-mV", type=float, default=-70.0)
    parser.add_argument("--equilibration-ms", type=float, default=500.0)
    parser.add_argument("--step-duration-ms", type=float, default=1000.0)
    parser.add_argument("--rin-step-duration-ms", type=float, default=300.0)
    parser.add_argument("--measure-start-ms", type=float, default=200.0)
    parser.add_argument("--measure-stop-ms", type=float, default=300.0)
    parser.add_argument("--post-ms", type=float, default=100.0)
    parser.add_argument("--dt-ms", type=float, default=0.05)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    parser.add_argument("--refractory-ms", type=float, default=2.0)
    parser.add_argument("--apply-slice-pruning", action="store_true")
    parser.add_argument("--soma-entry-depth-um", type=float, default=25.0)
    parser.add_argument("--slice-axis", choices=["x", "z"], default="z")
    parser.add_argument("--truncate-radius-um", type=float, default=500.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_mechanisms(args.repo_root.resolve())
    profiles = wang_ion_profiles(args.target_c)
    currents = parse_floats(args.current_values_nA)

    profile_rows = []
    rin_rows = []
    fi_rows = []

    for profile in profiles:
        profile_rows.append(
            {
                "profile": profile.name,
                "ena_mV": profile.ena_mV if profile.ena_mV is not None else "dura",
                "ek_mV": profile.ek_mV if profile.ek_mV is not None else "dura",
                "eca_mV": profile.eca_mV if profile.eca_mV is not None else "dura",
                "note": profile.note,
            }
        )

        rin_args = argparse.Namespace(**vars(args))
        rin_args.step_duration_ms = args.rin_step_duration_ms
        trace0 = simulate_step(args.cell_params, 0.0, profile, rin_args)
        trace_i = simulate_step(args.cell_params, args.test_current_nA, profile, rin_args)
        rest_v = baseline_mean(trace0, 0.0, 20.0)
        steady_v = baseline_mean(trace_i, args.measure_start_ms, args.measure_stop_ms)
        rin_rows.append(
            {
                "label": args.label,
                "profile": profile.name,
                "test_current_nA": args.test_current_nA,
                "rest_v_mV": rest_v,
                "steady_v_mV": steady_v,
                "delta_v_mV": steady_v - rest_v,
                "input_resistance_MOhm": (steady_v - rest_v) / args.test_current_nA,
                "retained_area_um2": trace_i["retained_area_um2"],
                "retained_area_fraction": trace_i["retained_area_fraction"],
            }
        )

        for current in currents:
            trace = simulate_step(args.cell_params, current, profile, args)
            spikes = spike_peaks(trace["t"], trace["v"], args.spike_threshold_mV, args.refractory_ms)
            spike_times = [t for t, _ in spikes if t <= args.step_duration_ms]
            during = trace["t"] <= args.step_duration_ms
            fi_rows.append(
                {
                    "label": args.label,
                    "profile": profile.name,
                    "current_nA": current,
                    "rest_v_mV": float(trace["v"][0]),
                    "max_v_mV": float(np.max(trace["v"][during])),
                    "min_v_mV": float(np.min(trace["v"][during])),
                    "spike_count_1s": len(spike_times),
                    "firing_rate_Hz": len(spike_times) / (args.step_duration_ms / 1000.0),
                    "first_spike_latency_ms": spike_times[0] if spike_times else float("nan"),
                    "retained_area_um2": trace["retained_area_um2"],
                    "retained_area_fraction": trace["retained_area_fraction"],
                }
            )

    with (args.out_dir / "ion_profiles.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(profile_rows[0].keys()))
        writer.writeheader()
        writer.writerows(profile_rows)
    with (args.out_dir / "input_resistance.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rin_rows[0].keys()))
        writer.writeheader()
        writer.writerows(rin_rows)
    with (args.out_dir / "current_injection_fi.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fi_rows[0].keys()))
        writer.writeheader()
        writer.writerows(fi_rows)

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.0), constrained_layout=True)
    names = [row["profile"] for row in rin_rows]
    x = np.arange(len(names))
    axes[0].bar(x, [row["input_resistance_MOhm"] for row in rin_rows], color="#4c78a8")
    axes[0].axhline(126.0, color="#d62728", ls="--", lw=1, label="Wang ChR2+")
    axes[0].axhline(110.0, color="#777777", ls=":", lw=1, label="Wang ChR2-")
    axes[0].set_xticks(x, names, rotation=25, ha="right")
    axes[0].set_ylabel("Rin (MOhm)")
    axes[0].set_title("Input resistance")
    axes[0].legend(fontsize=8)

    axes[1].bar(x, [row["rest_v_mV"] for row in rin_rows], color="#f58518")
    axes[1].axhline(-58.3, color="#d62728", ls="--", lw=1)
    axes[1].set_xticks(x, names, rotation=25, ha="right")
    axes[1].set_ylabel("rest Vm (mV)")
    axes[1].set_title("Resting voltage")

    for profile in profiles:
        sub = [row for row in fi_rows if row["profile"] == profile.name]
        axes[2].plot([row["current_nA"] for row in sub], [row["spike_count_1s"] for row in sub], marker="o", label=profile.name)
    axes[2].axhline(24.0, color="#d62728", ls="--", lw=1, label="Wang ~24 Hz")
    axes[2].set_xlabel("current (nA)")
    axes[2].set_ylabel("spikes / 1 s")
    axes[2].set_title("Somatic f-I")
    axes[2].legend(fontsize=7)
    fig.savefig(args.out_dir / "wang_solution_reversal_sweep.png", dpi=220)
    plt.close(fig)

    for row in profile_rows:
        print(row)
    for row in rin_rows:
        print(row)


if __name__ == "__main__":
    main()
