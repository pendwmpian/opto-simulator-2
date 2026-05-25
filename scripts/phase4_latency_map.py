from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from neuron import h

from chr2_kinetics import ChR2FourState
from phase4_area_threshold import (
    OpticalParams,
    local_irradiance_from_shape,
    section_midpoint_xyz,
)
from run_chr2_single_cell_sweep import CellFromNetPyNE, section_kind


@dataclass(frozen=True)
class MapResponse:
    row: int
    col: int
    center_x_um: float
    center_z_um: float
    size_um: float
    spike_count: int
    first_spike_latency_ms: float | None
    peak_v_mV: float
    max_depol_mV: float
    peak_chr2_current_nA: float
    chr2_charge_nC: float


def source_points_square(size_um: float, center_x_um: float, center_z_um: float, pitch_um: float) -> np.ndarray:
    half = size_um / 2.0
    xs = np.arange(center_x_um - half + pitch_um / 2.0, center_x_um + half, pitch_um)
    zs = np.arange(center_z_um - half + pitch_um / 2.0, center_z_um + half, pitch_um)
    return np.asarray([(x, z) for z in zs for x in xs], dtype=float)


def soma_center_xz(cell: CellFromNetPyNE) -> tuple[float, float]:
    soma = cell.sections["soma"]
    pts = [
        (float(h.x3d(i, sec=soma)), float(h.z3d(i, sec=soma)))
        for i in range(int(h.n3d(sec=soma)))
    ]
    if not pts:
        return 0.0, 0.0
    return float(np.mean([p[0] for p in pts])), float(np.mean([p[1] for p in pts]))


def run_one(args, row: int, col: int, center_x: float, center_z: float) -> tuple[MapResponse, np.ndarray, np.ndarray]:
    h.load_file("stdrun.hoc")
    cell = CellFromNetPyNE(args.cell_params)
    optical = OpticalParams(
        surface_irradiance_mw_mm2=args.surface_irradiance_mw_mm2,
        mu_eff_mm_inv=args.mu_eff_mm_inv,
        lateral_sigma0_um=args.lateral_sigma0_um,
        lateral_spread_per_depth=args.lateral_spread_per_depth,
        source_sample_pitch_um=args.source_sample_pitch_um,
    )
    src = source_points_square(args.stim_size_um, center_x, center_z, args.source_sample_pitch_um)
    surface_y = cell.surface_y()

    clamps = []
    kinetic_entries = []
    total_current_nA = 0.0
    for sec_name, sec in cell.sections.items():
        kind = section_kind(sec_name)
        if kind == "axon":
            continue
        for seg in sec:
            x, y, z = section_midpoint_xyz(sec, seg.x)
            irr = local_irradiance_from_shape(x, y, z, surface_y, src, optical)
            area_cm2 = float(h.area(seg.x, sec=sec)) * 1.0e-8
            g_mS_cm2 = args.base_g_mS_cm2_per_mw * args.expression_scale * irr
            g_uS = g_mS_cm2 * area_cm2 * 1000.0
            amp_nA = g_uS * (args.e_chr2_mV - args.v_init_mV)
            if amp_nA <= 0:
                continue
            clamp = h.IClamp(sec(seg.x))
            clamp.delay = args.pulse_start_ms
            if args.chr2_model == "fixed":
                clamp.dur = args.pulse_dur_ms
                clamp.amp = amp_nA
            else:
                clamp.dur = 1e9
                clamp.amp = 0.0
                kinetic_entries.append(
                    {
                        "clamp": clamp,
                        "seg": sec(seg.x),
                        "model": ChR2FourState(),
                        "irradiance": irr,
                        "area_cm2": area_cm2,
                    }
                )
            clamps.append(clamp)
            total_current_nA += amp_nA

    soma = cell.sections["soma"]
    h.dt = args.dt_ms
    h.tstop = args.tstop_ms
    h.v_init = args.v_init_mV
    h.finitialize(args.v_init_mV)
    if args.chr2_model == "fixed":
        t_vec = h.Vector().record(h._ref_t)
        v_vec = h.Vector().record(soma(0.5)._ref_v)
        h.continuerun(args.tstop_ms)
        t = np.asarray(t_vec)
        v = np.asarray(v_vec)
        peak_current_nA = total_current_nA
        charge_nC = total_current_nA * args.pulse_dur_ms / 1000.0
    else:
        t_values = []
        v_values = []
        current_values = []
        while h.t <= args.tstop_ms + 1e-9:
            total_i = 0.0
            light_on = args.pulse_start_ms <= h.t < args.pulse_start_ms + args.pulse_dur_ms
            for entry in kinetic_entries:
                irr = entry["irradiance"] if light_on else 0.0
                model = entry["model"]
                model.step(args.dt_ms, irr)
                gbar = args.base_g_mS_cm2_per_mw * args.expression_scale
                density_uA_cm2 = model.current_density_uA_cm2(gbar, float(entry["seg"].v))
                amp_nA = density_uA_cm2 * entry["area_cm2"] * 1000.0
                entry["clamp"].amp = max(0.0, amp_nA)
                total_i += max(0.0, amp_nA)
            t_values.append(float(h.t))
            v_values.append(float(soma(0.5).v))
            current_values.append(total_i)
            h.fadvance()
        t = np.asarray(t_values)
        v = np.asarray(v_values)
        current = np.asarray(current_values)
        peak_current_nA = float(np.max(current))
        charge_nC = float(np.trapezoid(current, t) / 1000.0)
    after = t >= args.pulse_start_ms
    baseline_window = (t >= args.pulse_start_ms - 20.0) & (t < args.pulse_start_ms)
    baseline_v = float(np.mean(v[baseline_window])) if np.any(baseline_window) else args.v_init_mV
    spike_indices = np.flatnonzero((v[1:] >= args.spike_threshold_mV) & (v[:-1] < args.spike_threshold_mV)) + 1
    spike_times = t[spike_indices]
    spike_times = spike_times[spike_times >= args.pulse_start_ms]
    latency = None if len(spike_times) == 0 else float(spike_times[0] - args.pulse_start_ms)
    response = MapResponse(
        row=row,
        col=col,
        center_x_um=float(center_x),
        center_z_um=float(center_z),
        size_um=args.stim_size_um,
        spike_count=int(len(spike_times)),
        first_spike_latency_ms=latency,
        peak_v_mV=float(np.max(v[after])),
        max_depol_mV=float(np.max(v[after]) - baseline_v),
        peak_chr2_current_nA=float(peak_current_nA),
        chr2_charge_nC=float(charge_nC),
    )
    return response, t, v


def format_table(values: np.ndarray, missing: str = "fail") -> list[list[str]]:
    table = []
    for row in values:
        out = []
        for value in row:
            if np.isnan(value):
                out.append(missing)
            else:
                out.append(f"{value:.2f}")
        table.append(out)
    return table


def plot_heatmap(matrix: np.ndarray, title: str, cbar_label: str, out: Path, failure_mask: np.ndarray | None = None) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.8), constrained_layout=True)
    masked = np.ma.masked_invalid(matrix)
    im = ax.imshow(masked, origin="upper", cmap="viridis")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_xlabel("x grid index")
    ax.set_ylabel("z grid index")
    ax.set_title(title)
    for r in range(matrix.shape[0]):
        for c in range(matrix.shape[1]):
            label = "fail" if np.isnan(matrix[r, c]) else f"{matrix[r, c]:.1f}"
            ax.text(c, r, label, ha="center", va="center", color="white", fontsize=8)
    if failure_mask is not None:
        for r, c in zip(*np.where(failure_mask)):
            ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False, edgecolor="red", lw=2))
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_waveforms(
    traces: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]],
    latency: np.ndarray,
    out: Path,
    pulse_start_ms: float,
    pulse_dur_ms: float,
) -> None:
    n = latency.shape[0]
    fig, axes = plt.subplots(n, n, figsize=(11, 9), sharex=True, sharey=True, constrained_layout=True)
    for r in range(n):
        for c in range(n):
            ax = axes[r, c]
            t, v = traces[(r, c)]
            ax.plot(t, v, color="#222222", lw=0.9)
            ax.axvspan(pulse_start_ms, pulse_start_ms + pulse_dur_ms, color="#f2c14e", alpha=0.25, lw=0)
            if np.isnan(latency[r, c]):
                ax.set_title("fail", fontsize=8)
            else:
                ax.set_title(f"{latency[r, c]:.1f} ms", fontsize=8)
            ax.set_ylim(-80, 45)
            ax.tick_params(labelsize=7)
    fig.supxlabel("time (ms)")
    fig.supylabel("soma Vm (mV)")
    fig.savefig(out, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/phase4_latency_map_100um_step50"))
    parser.add_argument("--field-size-um", type=float, default=200.0)
    parser.add_argument("--grid-n", type=int, default=4)
    parser.add_argument("--stim-size-um", type=float, default=100.0)
    parser.add_argument("--surface-irradiance-mw-mm2", type=float, default=1.0)
    parser.add_argument("--expression-scale", type=float, default=1.0)
    parser.add_argument("--mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--lateral-sigma0-um", type=float, default=8.0)
    parser.add_argument("--lateral-spread-per-depth", type=float, default=0.10)
    parser.add_argument("--source-sample-pitch-um", type=float, default=2.0)
    parser.add_argument("--base-g-mS-cm2-per-mw", type=float, default=0.08)
    parser.add_argument("--chr2-model", choices=["fixed", "kinetic"], default="fixed")
    parser.add_argument("--pulse-start-ms", type=float, default=100.0)
    parser.add_argument("--pulse-dur-ms", type=float, default=20.0)
    parser.add_argument("--tstop-ms", type=float, default=250.0)
    parser.add_argument("--dt-ms", type=float, default=0.025)
    parser.add_argument("--v-init-mV", type=float, default=-70.0)
    parser.add_argument("--e-chr2-mV", type=float, default=0.0)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    args = parser.parse_args()

    h.load_file("stdrun.hoc")
    ref_cell = CellFromNetPyNE(args.cell_params)
    soma_x, soma_z = soma_center_xz(ref_cell)
    del ref_cell

    step = args.field_size_um / args.grid_n
    offsets = np.linspace(-args.field_size_um / 2 + step / 2, args.field_size_um / 2 - step / 2, args.grid_n)

    responses: list[MapResponse] = []
    traces: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
    for r, z_off in enumerate(offsets):
        for c, x_off in enumerate(offsets):
            response, t, v = run_one(args, r, c, soma_x + x_off, soma_z + z_off)
            responses.append(response)
            traces[(r, c)] = (t, v)
            print(
                f"row={r} col={c} center=({response.center_x_um:.1f},{response.center_z_um:.1f}) "
                f"spikes={response.spike_count} latency={response.first_spike_latency_ms} peak={response.peak_v_mV:.2f}"
            )

    n = args.grid_n
    latency = np.full((n, n), np.nan)
    peak = np.full((n, n), np.nan)
    current = np.full((n, n), np.nan)
    for res in responses:
        if res.first_spike_latency_ms is not None:
            latency[res.row, res.col] = res.first_spike_latency_ms
        peak[res.row, res.col] = res.peak_v_mV
        current[res.row, res.col] = res.peak_chr2_current_nA

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "latency_map_results.json").open("w") as f:
        json.dump(
            {
                "args": vars(args) | {"cell_params": str(args.cell_params), "out_dir": str(args.out_dir)},
                "tables": {
                    "latency_ms": format_table(latency),
                    "peak_v_mV": format_table(peak),
                    "peak_chr2_current_nA": format_table(current),
                },
                "responses": [asdict(r) for r in responses],
            },
            f,
            indent=2,
        )
    with (args.out_dir / "latency_map_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(responses[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in responses)

    plot_heatmap(latency, "First spike latency", "latency (ms)", args.out_dir / "latency_table.png", np.isnan(latency))
    plot_heatmap(peak, "Peak soma Vm", "peak Vm (mV)", args.out_dir / "peak_vm_table.png")
    plot_heatmap(current, "Peak ChR2 current", "current (nA)", args.out_dir / "chr2_current_table.png")
    plot_waveforms(traces, latency, args.out_dir / "waveform_table.png", args.pulse_start_ms, args.pulse_dur_ms)

    print(args.out_dir / "latency_table.png")
    print(args.out_dir / "peak_vm_table.png")
    print(args.out_dir / "waveform_table.png")


if __name__ == "__main__":
    main()
