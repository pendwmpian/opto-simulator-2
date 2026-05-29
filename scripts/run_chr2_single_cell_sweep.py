from __future__ import annotations

import argparse
import json
import math
import os
import pickle
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from neuron import h


@dataclass(frozen=True)
class SweepResult:
    irradiance_mw_mm2: float
    expression_scale: float
    target: str
    peak_v_mV: float
    max_depol_mV: float
    spike_count: int
    first_spike_latency_ms: float | None
    total_chr2_charge_nC: float


def load_pickle(path: Path):
    with path.open("rb") as f:
        try:
            return pickle.load(f)
        except UnicodeDecodeError:
            f.seek(0)
            return pickle.load(f, encoding="latin1")


def section_kind(name: str) -> str:
    if name.startswith("soma"):
        return "soma"
    if name.startswith("axon"):
        return "axon"
    if name.startswith("apic") or name.startswith("Adend"):
        return "apic"
    if name.startswith("dend") or name.startswith("Bdend"):
        return "basal"
    return "other"


class CellFromNetPyNE:
    def __init__(self, cell_params: Path, use_original_biophysics: bool = False):
        self.rule = load_pickle(cell_params)
        self.sections: dict[str, h.Section] = {}
        self.use_original_biophysics = use_original_biophysics
        self._create_sections()
        self._connect_sections()
        self._apply_biophysics()

    def _create_sections(self) -> None:
        h.pt3dconst(1)
        for name, sec_rule in self.rule["secs"].items():
            sec = h.Section(name=name)
            self.sections[name] = sec
            geom = sec_rule.get("geom", {})
            sec.L = float(geom.get("L", 20.0))
            sec.diam = float(geom.get("diam", 1.0))
            sec.Ra = float(geom.get("Ra", 150.0))
            sec.cm = float(geom.get("cm", 1.0))
            sec.nseg = max(1, int(geom.get("nseg", 1)))
            if "pt3d" in geom:
                h.pt3dclear(sec=sec)
                for x, y, z, diam in geom["pt3d"]:
                    h.pt3dadd(float(x), float(y), float(z), float(diam), sec=sec)

    def _connect_sections(self) -> None:
        for name, sec_rule in self.rule["secs"].items():
            topol = sec_rule.get("topol", {})
            if not topol:
                continue
            parent = self.sections[str(topol["parentSec"])]
            self.sections[name].connect(parent(float(topol.get("parentX", 1.0))), float(topol.get("childX", 0.0)))

    def _apply_biophysics(self) -> None:
        if self.use_original_biophysics:
            self._apply_original_biophysics()
            return
        for name, sec in self.sections.items():
            kind = section_kind(name)
            sec.insert("pas")
            sec.g_pas = 3.0e-5
            sec.e_pas = -70.0
            if kind in {"soma", "axon"}:
                sec.insert("hh")
                if kind == "soma":
                    sec.gnabar_hh = 0.12
                    sec.gkbar_hh = 0.036
                    sec.gl_hh = 3.0e-5
                    sec.el_hh = -70.0
                else:
                    sec.gnabar_hh = 0.30
                    sec.gkbar_hh = 0.08
                    sec.gl_hh = 3.0e-5
                    sec.el_hh = -70.0
            elif kind in {"apic", "basal"}:
                sec.insert("hh")
                sec.gnabar_hh = 0.015
                sec.gkbar_hh = 0.006
                sec.gl_hh = 2.0e-5
                sec.el_hh = -70.0

    def _apply_original_biophysics(self) -> None:
        for name, sec in self.sections.items():
            sec_rule = self.rule["secs"][name]
            for mech_name, params in sec_rule.get("mechs", {}).items():
                sec.insert(mech_name)
                for param_name, value in params.items():
                    self._set_mech_param(sec, mech_name, param_name, value)
            self._apply_ion_params(sec, sec_rule.get("ions", {}))

    @staticmethod
    def _apply_ion_params(sec, ions: dict) -> None:
        ion_attr = {
            "na": {"i": "nai", "o": "nao", "e": "ena"},
            "k": {"i": "ki", "o": "ko", "e": "ek"},
            "ca": {"i": "cai", "o": "cao", "e": "eca"},
        }
        for ion_name, params in ions.items():
            for param_name, attr in ion_attr.get(ion_name, {}).items():
                if param_name in params:
                    setattr(sec, attr, float(params[param_name]))

    @staticmethod
    def _set_mech_param(sec, mech_name: str, param_name: str, value) -> None:
        attr = f"{param_name}_{mech_name}"
        if isinstance(value, list):
            if len(value) == 0:
                return
            segments = list(sec)
            if len(value) == len(segments):
                for seg, item in zip(segments, value):
                    setattr(getattr(seg, mech_name), param_name, float(item))
            else:
                xs = np.linspace(0.0, 1.0, len(value))
                for seg in segments:
                    item = float(np.interp(float(seg.x), xs, value))
                    setattr(getattr(seg, mech_name), param_name, item)
            return
        try:
            setattr(sec, attr, float(value))
        except Exception:
            for seg in sec:
                setattr(getattr(seg, mech_name), param_name, float(value))

    def surface_y(self) -> float:
        ys = []
        for sec_rule in self.rule["secs"].values():
            for point in sec_rule.get("geom", {}).get("pt3d", []):
                ys.append(float(point[1]))
        return max(ys)

    def segment_records(self, target: str, attenuation_um: float) -> list[dict]:
        surface_y = self.surface_y()
        rows = []
        for name, sec in self.sections.items():
            kind = section_kind(name)
            if target == "soma" and kind != "soma":
                continue
            if target == "dendrite" and kind not in {"apic", "basal"}:
                continue
            if target == "apical" and kind != "apic":
                continue
            if target == "basal" and kind != "basal":
                continue
            if target == "all" and kind == "axon":
                continue
            if kind == "axon":
                continue
            for seg in sec:
                y = self._segment_y(sec, seg.x)
                depth = max(0.0, surface_y - y)
                atten = math.exp(-depth / attenuation_um)
                area_cm2 = float(h.area(seg.x, sec=sec)) * 1.0e-8
                rows.append({"sec": sec, "x": seg.x, "kind": kind, "atten": atten, "area_cm2": area_cm2})
        return rows

    @staticmethod
    def _segment_y(sec, x: float) -> float:
        n3d = int(h.n3d(sec=sec))
        if n3d == 0:
            return 0.0
        if n3d == 1:
            return float(h.y3d(0, sec=sec))
        target_arc = float(x) * float(sec.L)
        prev_arc = float(h.arc3d(0, sec=sec))
        prev_y = float(h.y3d(0, sec=sec))
        for i in range(1, n3d):
            arc = float(h.arc3d(i, sec=sec))
            y = float(h.y3d(i, sec=sec))
            if arc >= target_arc:
                denom = max(1.0e-9, arc - prev_arc)
                frac = min(1.0, max(0.0, (target_arc - prev_arc) / denom))
                return prev_y + frac * (y - prev_y)
            prev_arc = arc
            prev_y = y
        return prev_y


def add_chr2_current(
    cell: CellFromNetPyNE,
    target: str,
    irradiance_mw_mm2: float,
    expression_scale: float,
    pulse_start_ms: float,
    pulse_dur_ms: float,
    attenuation_um: float,
    base_g_mS_cm2_per_mw: float,
    e_chr2_mV: float,
    v_ref_mV: float,
) -> tuple[list, float]:
    clamps = []
    total_i_nA = 0.0
    for row in cell.segment_records(target, attenuation_um):
        local_irr = irradiance_mw_mm2 * row["atten"]
        g_mS_cm2 = base_g_mS_cm2_per_mw * expression_scale * local_irr
        g_uS = g_mS_cm2 * row["area_cm2"] * 1000.0
        amp_nA = g_uS * (e_chr2_mV - v_ref_mV)
        if amp_nA <= 0:
            continue
        clamp = h.IClamp(row["sec"](row["x"]))
        clamp.delay = pulse_start_ms
        clamp.dur = pulse_dur_ms
        clamp.amp = amp_nA
        clamps.append(clamp)
        total_i_nA += amp_nA
    return clamps, total_i_nA * pulse_dur_ms / 1000.0


def run_condition(args, target: str, irradiance: float, expression_scale: float):
    h.load_file("stdrun.hoc")
    cell = CellFromNetPyNE(args.cell_params)
    clamps, charge_nC = add_chr2_current(
        cell=cell,
        target=target,
        irradiance_mw_mm2=irradiance,
        expression_scale=expression_scale,
        pulse_start_ms=args.pulse_start_ms,
        pulse_dur_ms=args.pulse_dur_ms,
        attenuation_um=args.attenuation_um,
        base_g_mS_cm2_per_mw=args.base_g_mS_cm2_per_mw,
        e_chr2_mV=args.e_chr2_mV,
        v_ref_mV=args.v_init_mV,
    )

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
    result = SweepResult(
        irradiance_mw_mm2=irradiance,
        expression_scale=expression_scale,
        target=target,
        peak_v_mV=float(np.max(v[after])),
        max_depol_mV=float(np.max(v[after]) - baseline_v),
        spike_count=int(len(spike_times)),
        first_spike_latency_ms=latency,
        total_chr2_charge_nC=float(charge_nC),
    )
    return result, t, v


def plot_summary(results: list[SweepResult], traces: dict[str, tuple[np.ndarray, np.ndarray]], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    targets = sorted({r.target for r in results})
    exprs = sorted({r.expression_scale for r in results})

    for expr in exprs:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
        for target in targets:
            rows = sorted([r for r in results if r.expression_scale == expr and r.target == target], key=lambda r: r.irradiance_mw_mm2)
            xs = [r.irradiance_mw_mm2 for r in rows]
            depol = [r.max_depol_mV for r in rows]
            latency = [np.nan if r.first_spike_latency_ms is None else r.first_spike_latency_ms for r in rows]
            axes[0].plot(xs, depol, marker="o", label=target)
            axes[1].plot(xs, latency, marker="o", label=target)
        axes[0].set_xscale("log")
        axes[0].set_xlabel("surface irradiance at 488 nm (mW/mm2)")
        axes[0].set_ylabel("max soma depolarization (mV)")
        axes[1].set_xscale("log")
        axes[1].set_xlabel("surface irradiance at 488 nm (mW/mm2)")
        axes[1].set_ylabel("first spike latency (ms)")
        axes[1].set_ylim(bottom=0)
        axes[0].legend(frameon=False)
        fig.suptitle(f"ChR2 expression scale = {expr:g}")
        fig.savefig(out_dir / f"sweep_expr_{expr:g}.png", dpi=200)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for label, (t, v) in traces.items():
        ax.plot(t, v, lw=1.2, label=label)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("soma Vm (mV)")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.savefig(out_dir / "example_traces.png", dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-params", type=Path, default=Path("external/M1_NetPyNE_CellReports_2023/sim/cells/PT5B_full_cellParams.pkl"))
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/phase1_chr2"))
    parser.add_argument("--irradiance", nargs="+", type=float, default=[0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
    parser.add_argument("--expression-scales", nargs="+", type=float, default=[0.25, 1.0, 4.0])
    parser.add_argument("--targets", nargs="+", default=["soma", "dendrite", "apical", "basal", "all"])
    parser.add_argument("--attenuation-um", type=float, default=471.7)
    parser.add_argument("--base-g-mS-cm2-per-mw", type=float, default=0.08)
    parser.add_argument("--pulse-start-ms", type=float, default=100.0)
    parser.add_argument("--pulse-dur-ms", type=float, default=20.0)
    parser.add_argument("--tstop-ms", type=float, default=250.0)
    parser.add_argument("--dt-ms", type=float, default=0.025)
    parser.add_argument("--v-init-mV", type=float, default=-70.0)
    parser.add_argument("--e-chr2-mV", type=float, default=0.0)
    parser.add_argument("--spike-threshold-mV", type=float, default=0.0)
    args = parser.parse_args()

    results: list[SweepResult] = []
    traces: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for expr in args.expression_scales:
        for target in args.targets:
            for irr in args.irradiance:
                result, t, v = run_condition(args, target, irr, expr)
                results.append(result)
                if expr in {args.expression_scales[0], args.expression_scales[-1]} and irr in {args.irradiance[0], args.irradiance[-1]} and target in {"soma", "dendrite", "all"}:
                    traces[f"{target}, {irr:g} mW/mm2, expr {expr:g}"] = (t, v)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "sweep_results.json").open("w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)
    plot_summary(results, traces, args.out_dir)
    for r in results:
        latency = "None" if r.first_spike_latency_ms is None else f"{r.first_spike_latency_ms:.2f}"
        print(
            f"{r.target:8s} expr={r.expression_scale:g} irr={r.irradiance_mw_mm2:g} "
            f"peak={r.peak_v_mV:.2f} depol={r.max_depol_mV:.2f} "
            f"spikes={r.spike_count} latency_ms={latency} charge_nC={r.total_chr2_charge_nC:.4f}"
        )


if __name__ == "__main__":
    main()
