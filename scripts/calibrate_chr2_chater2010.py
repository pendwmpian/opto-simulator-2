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


@dataclass(frozen=True)
class FoutzNikolicRates:
    """Foutz/Nikolic-style four-state rates in 1/ms.

    Values are converted from the commonly cited Foutz/Nikolic ChR2 parameter
    set where Kd1=130/s, Kd2=25/s, e12(light)=53/s, e21(light)=23/s,
    e12(dark)=22/s, e21(dark)=11/s, Kr=0.4/s.

    The light opening rates are retained from the previous simulator scaffold
    until the full photon-flux formulation is ported.
    """

    k1_per_mw_ms: float = 0.35
    k2_per_mw_ms: float = 0.12
    gd1_ms: float = 0.130
    gd2_ms: float = 0.025
    e12_light_ms: float = 0.053
    e21_light_ms: float = 0.023
    e12_dark_ms: float = 0.022
    e21_dark_ms: float = 0.011
    gr_ms: float = 0.0004
    gamma: float = 0.05

    def temperature_scaled(self, q10: float, reference_c: float = 37.0, target_c: float = 22.0) -> "FoutzNikolicRates":
        factor = q10 ** ((reference_c - target_c) / 10.0)
        return self.temperature_scaled_groups(q10, q10, q10, q10, reference_c, target_c)

    def temperature_scaled_groups(
        self,
        q10_open: float,
        q10_close: float,
        q10_adapt: float,
        q10_recovery: float,
        reference_c: float = 37.0,
        target_c: float = 22.0,
    ) -> "FoutzNikolicRates":
        exponent = (reference_c - target_c) / 10.0
        open_factor = q10_open**exponent
        close_factor = q10_close**exponent
        adapt_factor = q10_adapt**exponent
        recovery_factor = q10_recovery**exponent
        return FoutzNikolicRates(
            k1_per_mw_ms=self.k1_per_mw_ms / open_factor,
            k2_per_mw_ms=self.k2_per_mw_ms / open_factor,
            gd1_ms=self.gd1_ms / close_factor,
            gd2_ms=self.gd2_ms / close_factor,
            e12_light_ms=self.e12_light_ms / adapt_factor,
            e21_light_ms=self.e21_light_ms / adapt_factor,
            e12_dark_ms=self.e12_dark_ms / adapt_factor,
            e21_dark_ms=self.e21_dark_ms / adapt_factor,
            gr_ms=self.gr_ms / recovery_factor,
            gamma=self.gamma,
        )


@dataclass(frozen=True)
class ScaleSet:
    s_open: float
    s_close: float
    s_adapt: float
    s_recovery: float


@dataclass(frozen=True)
class ChaterMetrics:
    s_open: float
    s_close: float
    s_adapt: float
    s_recovery: float
    score: float
    flash_2ms_ttp_ms: float
    flash_2ms_decay_tau_ms: float | None
    flash_1ms_decay_tau_ms: float | None
    stim_300ms_decay_tau_ms: float | None
    stim_300ms_desensitization_fraction: float
    recovery_half_time_s: float | None
    q10: float | None = None
    q10_open: float | None = None
    q10_close: float | None = None
    q10_adapt: float | None = None
    q10_recovery: float | None = None


class FourState:
    def __init__(self, rates: FoutzNikolicRates, scales: ScaleSet):
        self.r = rates
        self.s = scales
        self.c1 = 1.0
        self.o1 = 0.0
        self.o2 = 0.0
        self.c2 = 0.0

    def open_fraction(self) -> float:
        return self.o1 + self.r.gamma * self.o2

    def step(self, dt_ms: float, irradiance_mw_mm2: float) -> None:
        light_on = irradiance_mw_mm2 > 0
        k1 = self.r.k1_per_mw_ms * self.s.s_open * irradiance_mw_mm2
        k2 = self.r.k2_per_mw_ms * self.s.s_open * irradiance_mw_mm2
        gd1 = self.r.gd1_ms * self.s.s_close
        gd2 = self.r.gd2_ms * self.s.s_close
        if light_on:
            e12 = self.r.e12_light_ms * self.s.s_adapt
            e21 = self.r.e21_light_ms * self.s.s_adapt
        else:
            e12 = self.r.e12_dark_ms * self.s.s_adapt
            e21 = self.r.e21_dark_ms * self.s.s_adapt
        gr = self.r.gr_ms * self.s.s_recovery

        dc1 = -k1 * self.c1 + gd1 * self.o1 + gr * self.c2
        do1 = k1 * self.c1 - gd1 * self.o1 - e12 * self.o1 + e21 * self.o2
        do2 = k2 * self.c2 + e12 * self.o1 - e21 * self.o2 - gd2 * self.o2
        dc2 = gd2 * self.o2 - k2 * self.c2 - gr * self.c2

        self.c1 = max(0.0, self.c1 + dt_ms * dc1)
        self.o1 = max(0.0, self.o1 + dt_ms * do1)
        self.o2 = max(0.0, self.o2 + dt_ms * do2)
        self.c2 = max(0.0, self.c2 + dt_ms * dc2)
        total = self.c1 + self.o1 + self.o2 + self.c2
        if total > 0:
            self.c1 /= total
            self.o1 /= total
            self.o2 /= total
            self.c2 /= total


def simulate_pulse(scales: ScaleSet, irradiance: float, pulse_ms: float, tstop_ms: float, dt_ms: float, rates: FoutzNikolicRates):
    model = FourState(rates, scales)
    times = np.arange(0.0, tstop_ms + dt_ms, dt_ms)
    current = np.zeros_like(times)
    for i, t in enumerate(times):
        irr = irradiance if 0.0 <= t < pulse_ms else 0.0
        model.step(dt_ms, irr)
        current[i] = model.open_fraction()
    if np.max(current) > 0:
        current = current / np.max(current)
    return times, current


def fit_decay_tau(t, y, start_ms, end_ms):
    mask = (t >= start_ms) & (t <= end_ms)
    if np.count_nonzero(mask) < 8:
        return None
    tt = t[mask] - start_ms
    yy = y[mask]
    baseline = float(np.min(yy[-max(3, len(yy) // 5) :]))
    amp = yy - baseline
    ok = amp > max(1e-6, 0.05 * np.max(amp))
    if np.count_nonzero(ok) < 8:
        return None
    slope, _ = np.polyfit(tt[ok], np.log(amp[ok]), 1)
    if slope >= 0:
        return None
    return float(-1.0 / slope)


def paired_pulse_recovery_half_time(scales: ScaleSet, irradiance: float, dt_ms: float, rates: FoutzNikolicRates):
    intervals_ms = [50, 100, 200, 500, 1000, 2000, 3100, 5000, 8000]
    ratios = []
    for interval in intervals_ms:
        model = FourState(rates, scales)
        tstop = 200.0 + interval + 20.0
        times = np.arange(0.0, tstop + dt_ms, dt_ms)
        first_peak = 0.0
        second_peak = 0.0
        for t in times:
            light = (0.0 <= t < 200.0) or (200.0 + interval <= t < 202.0 + interval)
            model.step(dt_ms, irradiance if light else 0.0)
            cur = model.open_fraction()
            if 0.0 <= t < 200.0:
                first_peak = max(first_peak, cur)
            if 200.0 + interval <= t < 220.0 + interval:
                second_peak = max(second_peak, cur)
        ratios.append(second_peak / first_peak if first_peak > 0 else 0.0)
    target = 0.5
    for i in range(1, len(intervals_ms)):
        if ratios[i - 1] <= target <= ratios[i]:
            x0, x1 = intervals_ms[i - 1], intervals_ms[i]
            y0, y1 = ratios[i - 1], ratios[i]
            if y1 == y0:
                return float(x1 / 1000.0)
            frac = (target - y0) / (y1 - y0)
            return float((x0 + frac * (x1 - x0)) / 1000.0)
    return None


def evaluate(scales: ScaleSet, rates: FoutzNikolicRates, irradiance: float, dt_ms: float) -> ChaterMetrics:
    t2, y2 = simulate_pulse(scales, irradiance, 2.0, 80.0, dt_ms, rates)
    ttp2 = float(t2[int(np.argmax(y2))])
    tau2 = fit_decay_tau(t2, y2, 2.0, 40.0)

    t1, y1 = simulate_pulse(scales, irradiance, 1.0, 80.0, dt_ms, rates)
    tau1 = fit_decay_tau(t1, y1, 1.0, 40.0)

    t300, y300 = simulate_pulse(scales, irradiance, 300.0, 450.0, dt_ms, rates)
    tau300 = fit_decay_tau(t300, y300, 300.0, 380.0)
    peak_300 = float(np.max(y300[t300 < 300.0]))
    end_300 = float(y300[np.where(t300 < 300.0)[0][-1]])
    desens = 1.0 - end_300 / peak_300 if peak_300 > 0 else 0.0

    rec = paired_pulse_recovery_half_time(scales, irradiance, dt_ms, rates)

    targets = {
        "ttp2": (2.9, 0.3),
        "tau2": (7.2, 1.0),
        "tau1": (10.7, 1.0),
        "tau300": (11.5, 1.5),
        "desens": (0.75, 0.15),
        "rec": (3.1, 0.8),
    }
    score = ((ttp2 - targets["ttp2"][0]) / targets["ttp2"][1]) ** 2
    if tau2 is not None:
        score += ((tau2 - targets["tau2"][0]) / targets["tau2"][1]) ** 2
    else:
        score += 100.0
    if tau1 is not None:
        score += ((tau1 - targets["tau1"][0]) / targets["tau1"][1]) ** 2
    else:
        score += 100.0
    if tau300 is not None:
        score += ((tau300 - targets["tau300"][0]) / targets["tau300"][1]) ** 2
    else:
        score += 100.0
    score += ((desens - targets["desens"][0]) / targets["desens"][1]) ** 2
    if rec is not None:
        score += ((rec - targets["rec"][0]) / targets["rec"][1]) ** 2
    else:
        score += 100.0

    return ChaterMetrics(
        s_open=scales.s_open,
        s_close=scales.s_close,
        s_adapt=scales.s_adapt,
        s_recovery=scales.s_recovery,
        score=float(score),
        flash_2ms_ttp_ms=ttp2,
        flash_2ms_decay_tau_ms=tau2,
        flash_1ms_decay_tau_ms=tau1,
        stim_300ms_decay_tau_ms=tau300,
        stim_300ms_desensitization_fraction=desens,
        recovery_half_time_s=rec,
    )


def evaluate_q10(q10: float, rates: FoutzNikolicRates, irradiance: float, dt_ms: float) -> ChaterMetrics:
    scaled_rates = rates.temperature_scaled(q10)
    metrics = evaluate(ScaleSet(1.0, 1.0, 1.0, 1.0), scaled_rates, irradiance, dt_ms)
    return ChaterMetrics(**(asdict(metrics) | {"q10": q10}))


def evaluate_q10_groups(
    q10_open: float,
    q10_close: float,
    q10_adapt: float,
    q10_recovery: float,
    rates: FoutzNikolicRates,
    irradiance: float,
    dt_ms: float,
) -> ChaterMetrics:
    scaled_rates = rates.temperature_scaled_groups(q10_open, q10_close, q10_adapt, q10_recovery)
    metrics = evaluate(ScaleSet(1.0, 1.0, 1.0, 1.0), scaled_rates, irradiance, dt_ms)
    return ChaterMetrics(
        **(
            asdict(metrics)
            | {
                "q10_open": q10_open,
                "q10_close": q10_close,
                "q10_adapt": q10_adapt,
                "q10_recovery": q10_recovery,
            }
        )
    )


def plot_best(best: ChaterMetrics, rates: FoutzNikolicRates, irradiance: float, dt_ms: float, out: Path):
    scales = ScaleSet(best.s_open, best.s_close, best.s_adapt, best.s_recovery)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), constrained_layout=True)
    for ax, pulse, title in [(axes[0], 1.0, "1 ms flash"), (axes[1], 2.0, "2 ms flash"), (axes[2], 300.0, "300 ms stimulus")]:
        tstop = 80.0 if pulse < 10 else 450.0
        t, y = simulate_pulse(scales, irradiance, pulse, tstop, dt_ms, rates)
        ax.plot(t, y, color="#222222", lw=1.3)
        ax.axvspan(0, pulse, color="#f0c419", alpha=0.25)
        ax.set_xlabel("time (ms)")
        ax.set_ylabel("normalized current")
        ax.set_title(title)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/chater2010_chr2_temperature"))
    parser.add_argument("--irradiance-mw-mm2", type=float, default=1.0)
    parser.add_argument("--dt-ms", type=float, default=0.02)
    parser.add_argument("--mode", choices=["grid", "q10", "q10-groups"], default="grid")
    parser.add_argument("--q10-values", default="1.5,2.0,2.5,3.0")
    args = parser.parse_args()

    rates = FoutzNikolicRates()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "q10":
        results = [
            evaluate_q10(float(raw.strip()), rates, args.irradiance_mw_mm2, args.dt_ms)
            for raw in args.q10_values.split(",")
            if raw.strip()
        ]
        results.sort(key=lambda r: r.score)
        best = results[0]
        with (args.out_dir / "chater2010_q10_results.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(asdict(best).keys()))
            writer.writeheader()
            writer.writerows(asdict(r) for r in results)
        with (args.out_dir / "chater2010_q10_best.json").open("w") as f:
            json.dump({"baseline_rates": asdict(rates), "best": asdict(best)}, f, indent=2)
        best_rates = rates.temperature_scaled(best.q10 if best.q10 is not None else 2.0)
        plot_best(best, best_rates, args.irradiance_mw_mm2, args.dt_ms, args.out_dir / "chater2010_q10_best_traces.png")
        print(best)
        return
    if args.mode == "q10-groups":
        values = [float(raw.strip()) for raw in args.q10_values.split(",") if raw.strip()]
        results = []
        for q_open in values:
            for q_close in values:
                for q_adapt in values:
                    for q_recovery in values:
                        results.append(
                            evaluate_q10_groups(
                                q_open,
                                q_close,
                                q_adapt,
                                q_recovery,
                                rates,
                                args.irradiance_mw_mm2,
                                args.dt_ms,
                            )
                        )
        results.sort(key=lambda r: r.score)
        best = results[0]
        with (args.out_dir / "chater2010_q10_group_results.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(asdict(best).keys()))
            writer.writeheader()
            writer.writerows(asdict(r) for r in results)
        with (args.out_dir / "chater2010_q10_group_best.json").open("w") as f:
            json.dump({"baseline_rates": asdict(rates), "best": asdict(best)}, f, indent=2)
        best_rates = rates.temperature_scaled_groups(
            best.q10_open or 2.0,
            best.q10_close or 2.0,
            best.q10_adapt or 2.0,
            best.q10_recovery or 2.0,
        )
        plot_best(best, best_rates, args.irradiance_mw_mm2, args.dt_ms, args.out_dir / "chater2010_q10_group_best_traces.png")
        print(best)
        return

    scale_values = {
        "open": [0.25, 0.35, 0.5, 0.7, 1.0],
        "close": [0.5, 0.7, 1.0, 1.4, 2.0],
        "adapt": [0.5, 0.7, 1.0, 1.4, 2.0],
        "recovery": [0.1, 0.2, 0.4, 0.7, 1.0],
    }
    results = []
    for s_open in scale_values["open"]:
        for s_close in scale_values["close"]:
            for s_adapt in scale_values["adapt"]:
                for s_recovery in scale_values["recovery"]:
                    scales = ScaleSet(s_open, s_close, s_adapt, s_recovery)
                    results.append(evaluate(scales, rates, args.irradiance_mw_mm2, args.dt_ms))
    results.sort(key=lambda r: r.score)
    best = results[0]

    with (args.out_dir / "chater2010_grid_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(best).keys()))
        writer.writeheader()
        writer.writerows(asdict(r) for r in results)
    with (args.out_dir / "chater2010_best.json").open("w") as f:
        json.dump({"rates": asdict(rates), "best": asdict(best)}, f, indent=2)
    plot_best(best, rates, args.irradiance_mw_mm2, args.dt_ms, args.out_dir / "chater2010_best_traces.png")
    print(best)


if __name__ == "__main__":
    main()
