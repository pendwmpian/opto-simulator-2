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

from run_chr2_single_cell_sweep import CellFromNetPyNE, section_kind


@dataclass(frozen=True)
class OpticalParams:
    surface_irradiance_mw_mm2: float
    mu_eff_mm_inv: float
    lateral_sigma0_um: float
    lateral_spread_per_depth: float
    source_sample_pitch_um: float


@dataclass(frozen=True)
class Response:
    shape: str
    size_um: float
    irradiance_mw_mm2: float
    expression_scale: float
    spike_count: int
    first_spike_latency_ms: float | None
    peak_v_mV: float
    max_depol_mV: float
    peak_chr2_current_nA: float
    chr2_charge_nC: float


def section_midpoint_xyz(sec, x: float) -> tuple[float, float, float]:
    n3d = int(h.n3d(sec=sec))
    if n3d == 0:
        return 0.0, 0.0, 0.0
    if n3d == 1:
        return float(h.x3d(0, sec=sec)), float(h.y3d(0, sec=sec)), float(h.z3d(0, sec=sec))
    target_arc = float(x) * float(sec.L)
    prev_arc = float(h.arc3d(0, sec=sec))
    prev = np.array([h.x3d(0, sec=sec), h.y3d(0, sec=sec), h.z3d(0, sec=sec)], dtype=float)
    for i in range(1, n3d):
        arc = float(h.arc3d(i, sec=sec))
        cur = np.array([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)], dtype=float)
        if arc >= target_arc:
            frac = min(1.0, max(0.0, (target_arc - prev_arc) / max(1.0e-9, arc - prev_arc)))
            p = prev + frac * (cur - prev)
            return float(p[0]), float(p[1]), float(p[2])
        prev_arc = arc
        prev = cur
    return float(prev[0]), float(prev[1]), float(prev[2])


def soma_center_xz(cell: CellFromNetPyNE) -> tuple[float, float]:
    soma = cell.sections["soma"]
    pts = [
        (float(h.x3d(i, sec=soma)), float(h.z3d(i, sec=soma)))
        for i in range(int(h.n3d(sec=soma)))
    ]
    if not pts:
        return 0.0, 0.0
    return float(np.mean([p[0] for p in pts])), float(np.mean([p[1] for p in pts]))


def source_points(shape: str, size_um: float, center_x_um: float, center_z_um: float, pitch_um: float) -> np.ndarray:
    half = size_um / 2.0
    xs = np.arange(center_x_um - half + pitch_um / 2, center_x_um + half, pitch_um)
    zs = np.arange(center_z_um - half + pitch_um / 2, center_z_um + half, pitch_um)
    points = []
    radius = size_um / 2.0
    for z in zs:
        for x in xs:
            if shape == "square":
                points.append((x, z))
            elif shape == "circle":
                if (x - center_x_um) ** 2 + (z - center_z_um) ** 2 <= radius**2:
                    points.append((x, z))
            else:
                raise ValueError(shape)
    if not points:
        points = [(center_x_um, center_z_um)]
    return np.asarray(points, dtype=float)


def local_irradiance_from_shape(
    x_um: float,
    y_um: float,
    z_um: float,
    surface_y_um: float,
    src: np.ndarray,
    optical: OpticalParams,
) -> float:
    depth_um = max(0.0, surface_y_um - y_um)
    attenuation = math.exp(-optical.mu_eff_mm_inv * depth_um / 1000.0)
    sigma_um = max(1.0, optical.lateral_sigma0_um + optical.lateral_spread_per_depth * depth_um)
    pixel_area_um2 = optical.source_sample_pitch_um**2
    dx = x_um - src[:, 0]
    dz = z_um - src[:, 1]
    r2 = dx * dx + dz * dz
    weights = pixel_area_um2 / (2.0 * math.pi * sigma_um * sigma_um) * np.exp(-0.5 * r2 / (sigma_um * sigma_um))
    return float(optical.surface_irradiance_mw_mm2 * attenuation * np.sum(weights))


def segment_inputs(cell: CellFromNetPyNE, shape: str, size_um: float, optical: OpticalParams):
    center_x, center_z = soma_center_xz(cell)
    src = source_points(shape, size_um, center_x, center_z, optical.source_sample_pitch_um)
    surface_y = cell.surface_y()
    rows = []
    for sec_name, sec in cell.sections.items():
        kind = section_kind(sec_name)
        if kind == "axon":
            continue
        for seg in sec:
            x, y, z = section_midpoint_xyz(sec, seg.x)
            area_cm2 = float(h.area(seg.x, sec=sec)) * 1.0e-8
            irr = local_irradiance_from_shape(x, y, z, surface_y, src, optical)
            rows.append({"sec": sec, "x": seg.x, "kind": kind, "irradiance": irr, "area_cm2": area_cm2})
    return rows


def run_response(args, shape: str, size_um: float) -> tuple[Response, np.ndarray, np.ndarray]:
    h.load_file("stdrun.hoc")
    cell = CellFromNetPyNE(args.cell_params)
    optical = OpticalParams(
        surface_irradiance_mw_mm2=args.surface_irradiance_mw_mm2,
        mu_eff_mm_inv=args.mu_eff_mm_inv,
        lateral_sigma0_um=args.lateral_sigma0_um,
        lateral_spread_per_depth=args.lateral_spread_per_depth,
        source_sample_pitch_um=args.source_sample_pitch_um,
    )
    inputs = segment_inputs(cell, shape, size_um, optical)

    clamps = []
    total_current_nA = 0.0
    for row in inputs:
        g_mS_cm2 = args.base_g_mS_cm2_per_mw * args.expression_scale * row["irradiance"]
        g_uS = g_mS_cm2 * row["area_cm2"] * 1000.0
        amp_nA = g_uS * (args.e_chr2_mV - args.v_init_mV)
        if amp_nA <= 0:
            continue
        clamp = h.IClamp(row["sec"](row["x"]))
        clamp.delay = args.pulse_start_ms
        clamp.dur = args.pulse_dur_ms
        clamp.amp = amp_nA
        clamps.append(clamp)
        total_current_nA += amp_nA

    soma = cell.sections["soma"]
    t_vec = h.Vector().record(h._ref_t)
    v_vec = h.Vector().record(soma(0.5)._ref_v)

    h.dt = args.dt_ms
    h.tstop = args.tstop_ms
    h.v_init = args.v_init_mV
    h.finitialize(args.v_init_mV)
    h.continuerun(args.tstop_ms)

    t = np.asarray(t_vec)
    v = np.asarray(v_vec)
    after = t >= args.pulse_start_ms
    baseline_window = (t >= args.pulse_start_ms - 20.0) & (t < args.pulse_start_ms)
    baseline_v = float(np.mean(v[baseline_window])) if np.any(baseline_window) else args.v_init_mV
    spike_indices = np.flatnonzero((v[1:] >= args.spike_threshold_mV) & (v[:-1] < args.spike_threshold_mV)) + 1
    spike_times = t[spike_indices]
    spike_times = spike_times[spike_times >= args.pulse_start_ms]
    latency = None if len(spike_times) == 0 else float(spike_times[0] - args.pulse_start_ms)
    response = Response(
        shape=shape,
        size_um=size_um,
        irradiance_mw_mm2=args.surface_irradiance_mw_mm2,
        expression_scale=args.expression_scale,
        spike_count=int(len(spike_times)),
        first_spike_latency_ms=latency,
        peak_v_mV=float(np.max(v[after])),
        max_depol_mV=float(np.max(v[after]) - baseline_v),
        peak_chr2_current_nA=float(total_current_nA),
        chr2_charge_nC=float(total_current_nA * args.pulse_dur_ms / 1000.0),
    )
    return response, t, v


def binary_search_threshold(args, shape: str) -> tuple[float | None, list[Response], dict[str, tuple[np.ndarray, np.ndarray]]]:
    responses: list[Response] = []
    traces: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    lo = args.min_size_um
    hi = args.max_size_um
    r_lo, t_lo, v_lo = run_response(args, shape, lo)
    responses.append(r_lo)
    traces[f"{shape} {lo:g} um"] = (t_lo, v_lo)
    if r_lo.spike_count > 0:
        return lo, responses, traces

    r_hi, t_hi, v_hi = run_response(args, shape, hi)
    responses.append(r_hi)
    traces[f"{shape} {hi:g} um"] = (t_hi, v_hi)
    if r_hi.spike_count == 0:
        return None, responses, traces

    for _ in range(args.binary_iterations):
        mid = 0.5 * (lo + hi)
        r_mid, t_mid, v_mid = run_response(args, shape, mid)
        responses.append(r_mid)
        if r_mid.spike_count > 0:
            hi = mid
            traces[f"{shape} {mid:g} um"] = (t_mid, v_mid)
        else:
            lo = mid
    return hi, responses, traces


def plot_results(results: list[Response], traces: dict[str, tuple[np.ndarray, np.ndarray]], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    for shape in sorted({r.shape for r in results}):
        rows = sorted([r for r in results if r.shape == shape], key=lambda r: r.size_um)
        xs = [r.size_um for r in rows]
        axes[0].plot(xs, [r.max_depol_mV for r in rows], marker="o", label=shape)
        axes[1].plot(xs, [r.peak_chr2_current_nA for r in rows], marker="o", label=shape)
        axes[2].plot(
            xs,
            [np.nan if r.first_spike_latency_ms is None else r.first_spike_latency_ms for r in rows],
            marker="o",
            label=shape,
        )
    axes[0].set_xlabel("stimulus size (um)")
    axes[0].set_ylabel("max soma depolarization (mV)")
    axes[1].set_xlabel("stimulus size (um)")
    axes[1].set_ylabel("peak ChR2 photocurrent (nA)")
    axes[2].set_xlabel("stimulus size (um)")
    axes[2].set_ylabel("first spike latency (ms)")
    for ax in axes:
        ax.legend(frameon=False)
    fig.savefig(out_dir / "area_threshold_summary.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for label, (t, v) in traces.items():
        ax.plot(t, v, lw=1.1, label=label)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("soma Vm (mV)")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.savefig(out_dir / "area_threshold_traces.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/phase4_area_threshold"))
    parser.add_argument("--surface-irradiance-mw-mm2", type=float, default=1.0)
    parser.add_argument("--expression-scale", type=float, default=1.0)
    parser.add_argument("--mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--lateral-sigma0-um", type=float, default=8.0)
    parser.add_argument("--lateral-spread-per-depth", type=float, default=0.10)
    parser.add_argument("--source-sample-pitch-um", type=float, default=2.0)
    parser.add_argument("--min-size-um", type=float, default=2.0)
    parser.add_argument("--max-size-um", type=float, default=200.0)
    parser.add_argument("--binary-iterations", type=int, default=8)
    parser.add_argument("--base-g-mS-cm2-per-mw", type=float, default=0.08)
    parser.add_argument("--pulse-start-ms", type=float, default=100.0)
    parser.add_argument("--pulse-dur-ms", type=float, default=20.0)
    parser.add_argument("--tstop-ms", type=float, default=250.0)
    parser.add_argument("--dt-ms", type=float, default=0.025)
    parser.add_argument("--v-init-mV", type=float, default=-70.0)
    parser.add_argument("--e-chr2-mV", type=float, default=0.0)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    args = parser.parse_args()

    all_results: list[Response] = []
    all_traces: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    thresholds = {}
    for shape in ["square", "circle"]:
        threshold, responses, traces = binary_search_threshold(args, shape)
        thresholds[shape] = threshold
        all_results.extend(responses)
        all_traces.update(traces)
        if threshold is None:
            print(f"{shape}: no spike up to {args.max_size_um:g} um")
        else:
            print(f"{shape}: threshold ~= {threshold:.2f} um")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "area_threshold_results.json").open("w") as f:
        json.dump(
            {
                "thresholds_um": thresholds,
                "responses": [asdict(r) for r in all_results],
                "args": vars(args) | {"cell_params": str(args.cell_params), "out_dir": str(args.out_dir)},
            },
            f,
            indent=2,
        )
    with (args.out_dir / "area_threshold_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(all_results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in all_results)
    plot_results(all_results, all_traces, args.out_dir)
    print(args.out_dir / "area_threshold_summary.png")
    print(args.out_dir / "area_threshold_traces.png")


if __name__ == "__main__":
    main()
