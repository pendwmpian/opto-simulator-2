from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ChR2KineticsParams:
    """Compact four-state ChR2 kinetics.

    Units:
        irradiance_mw_mm2: mW/mm2
        rates: 1/ms
        gbar_mS_cm2_per_open: mS/cm2
        voltages: mV

    This is a first implementation scaffold following the Nikolic/Grossman/Foutz
    four-state structure. The rate constants are intentionally centralized so
    patch-clamp data can replace them without changing simulation code.
    """

    e_chr2_mV: float = 0.0
    gamma: float = 0.1
    k1_per_mw_ms: float = 0.35
    k2_per_mw_ms: float = 0.12
    gd1_ms: float = 0.10
    gd2_ms: float = 0.025
    e12_ms: float = 0.011
    e21_ms: float = 0.008
    gr_ms: float = 0.00033


@dataclass
class ChR2State:
    c1: float = 1.0
    o1: float = 0.0
    o2: float = 0.0
    c2: float = 0.0

    def normalize(self) -> None:
        total = self.c1 + self.o1 + self.o2 + self.c2
        if total <= 0:
            self.c1, self.o1, self.o2, self.c2 = 1.0, 0.0, 0.0, 0.0
            return
        self.c1 /= total
        self.o1 /= total
        self.o2 /= total
        self.c2 /= total


class ChR2FourState:
    """Four-state ChR2 photocycle with voltage-dependent driving force.

    State scheme:
        C1 --light--> O1 --desensitization--> O2 --close--> C2 --recover--> C1
        O1 and O2 can interconvert.

    The channel conductance factor is:
        open_fraction = O1 + gamma * O2

    The inward current density is:
        I_inward = gbar * open_fraction * irradiance-independent_state * (E - V)

    The caller converts current density to total segment current.
    """

    def __init__(self, params: ChR2KineticsParams | None = None):
        self.params = params or ChR2KineticsParams()
        self.state = ChR2State()

    def reset(self) -> None:
        self.state = ChR2State()

    def step(self, dt_ms: float, irradiance_mw_mm2: float) -> None:
        p = self.params
        s = self.state
        ka1 = p.k1_per_mw_ms * max(0.0, irradiance_mw_mm2)
        ka2 = p.k2_per_mw_ms * max(0.0, irradiance_mw_mm2)

        dc1 = -ka1 * s.c1 + p.gd1_ms * s.o1 + p.gr_ms * s.c2
        do1 = ka1 * s.c1 - p.gd1_ms * s.o1 - p.e12_ms * s.o1 + p.e21_ms * s.o2
        do2 = ka2 * s.c2 + p.e12_ms * s.o1 - p.e21_ms * s.o2 - p.gd2_ms * s.o2
        dc2 = p.gd2_ms * s.o2 - ka2 * s.c2 - p.gr_ms * s.c2

        s.c1 = max(0.0, s.c1 + dt_ms * dc1)
        s.o1 = max(0.0, s.o1 + dt_ms * do1)
        s.o2 = max(0.0, s.o2 + dt_ms * do2)
        s.c2 = max(0.0, s.c2 + dt_ms * dc2)
        s.normalize()

    def open_fraction(self) -> float:
        p = self.params
        s = self.state
        return s.o1 + p.gamma * s.o2

    def current_density_uA_cm2(self, gbar_mS_cm2: float, v_mV: float) -> float:
        """Return inward current density in uA/cm2."""
        driving_force_mV = self.params.e_chr2_mV - v_mV
        return gbar_mS_cm2 * self.open_fraction() * driving_force_mV

