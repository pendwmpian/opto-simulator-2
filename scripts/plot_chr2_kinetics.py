from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from chr2_kinetics import ChR2FourState


def simulate(dt_ms: float, tstop_ms: float, pulse_start_ms: float, pulse_dur_ms: float, irradiance: float, v_hold: float, gbar: float):
    model = ChR2FourState()
    t = np.arange(0.0, tstop_ms + dt_ms, dt_ms)
    open_fraction = np.zeros_like(t)
    current = np.zeros_like(t)
    light = np.zeros_like(t)
    for i, ti in enumerate(t):
        irr = irradiance if pulse_start_ms <= ti < pulse_start_ms + pulse_dur_ms else 0.0
        light[i] = irr
        model.step(dt_ms, irr)
        open_fraction[i] = model.open_fraction()
        current[i] = model.current_density_uA_cm2(gbar, v_hold)
    return t, light, open_fraction, current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("outputs/chr2_kinetics/chr2_50ms_voltage_clamp.png"))
    parser.add_argument("--dt-ms", type=float, default=0.025)
    parser.add_argument("--tstop-ms", type=float, default=250.0)
    parser.add_argument("--pulse-start-ms", type=float, default=50.0)
    parser.add_argument("--pulse-dur-ms", type=float, default=50.0)
    parser.add_argument("--gbar-mS-cm2", type=float, default=0.08)
    parser.add_argument("--v-hold-mV", type=float, default=-70.0)
    parser.add_argument("--irradiance", nargs="+", type=float, default=[0.1, 0.5, 1.0, 5.0])
    args = parser.parse_args()

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, constrained_layout=True)
    for irr in args.irradiance:
        t, light, open_fraction, current = simulate(
            args.dt_ms,
            args.tstop_ms,
            args.pulse_start_ms,
            args.pulse_dur_ms,
            irr,
            args.v_hold_mV,
            args.gbar_mS_cm2,
        )
        axes[0].plot(t, open_fraction, label=f"{irr:g} mW/mm2")
        axes[1].plot(t, current, label=f"{irr:g} mW/mm2")

    for ax in axes:
        ax.axvspan(args.pulse_start_ms, args.pulse_start_ms + args.pulse_dur_ms, color="#f2c14e", alpha=0.25, lw=0)
        ax.legend(frameon=False)
    axes[0].set_ylabel("open fraction")
    axes[1].set_ylabel("inward current density (uA/cm2)")
    axes[1].set_xlabel("time (ms)")
    fig.suptitle(f"Four-state ChR2 kinetics, Vhold={args.v_hold_mV:g} mV, gbar={args.gbar_mS_cm2:g} mS/cm2")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    plt.close(fig)
    print(args.out)


if __name__ == "__main__":
    main()
