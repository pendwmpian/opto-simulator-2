# PyRhO six-state Wang Fig. 2A-C implementation plan

## Goal

Replace the Foutz/Nikolic four-state ChR2 mechanism for Wang Fig. 2A-C calibration with a PyRhO/Grossman-style six-state ChR2 photocycle.

For this first step, only reproduce Wang Fig. 2A-C-style voltage-clamp behavior:

- estimate `gbar`,
- simulate 1-s photocurrents,
- simulate intensity-response curves,
- report peak, steady current, time-to-peak, inactivation tau, and intensity half-max.

Fig. 2D-F flash-duration simulations are intentionally deferred.

## Six-state model structure

Use a PyRhO/Grossman-style six-state model:

`C1 -> I1 -> O1 <-> O2 <- I2 <- C2`, with recovery/desensitization links.

State meanings:

- `C1`: dark-adapted closed state.
- `I1`: activation intermediate after photon absorption from `C1`.
- `O1`: high-conductance open state.
- `O2`: low-conductance/light-adapted open state.
- `I2`: activation intermediate from `C2`.
- `C2`: light-adapted/desensitized closed state.

Conducting occupancy:

`open = O1 + gamma * O2`

where `gamma` is the relative conductance of `O2`.

## Parameter assumptions

Use PyRhO six-state ChR2 fit parameters from Evans et al. 2016 / PyRhO:

- `phi_m = 5.02e17 photons/mm2/s`
- `k1 = 18.2 /ms`
- `k2 = 4.07 /ms`
- `p = 0.981`
- `Gf0 = 0.0365 /ms`
- `kf = 0.121 /ms`
- `Gb0 = 0.0143 /ms`
- `kb = 0.131 /ms`
- `q = 1.45`
- `Go1 = 1.93 /ms`
- `Go2 = 3.38 /ms`
- `Gd1 = 0.108 /ms`
- `Gd2 = 0.0115 /ms`
- `Gr0 = 0.00033 /ms`
- `gamma = 0.00369`
- `E_ChR2 = 0 mV`

`gbar` is not taken from PyRhO `g0`. Instead, `gbar` remains the membrane-area normalized maximum conductance in `mS/cm2`, so it can be compared with the previous model and calibrated to Wang line18 data.

## Light conversion

The NEURON mechanism receives irradiance in `mW/mm2`.

It converts to photon flux as:

`phi = 1000 * irradiance / E_photon / 1e6`

where `E_photon = hc / wavelength`, and the final unit is `photons/mm2/s`.

## Q10 sweep

Because the six-state model cannot use Williams parameter-specific Q10 mapping directly, use uniform kinetic scaling:

- `Q10 = 1.0`: no kinetic correction.
- `Q10 = 1.5`
- `Q10 = 2.0`

For Q10 > 1, all kinetic rates are scaled from a 37 deg C reference to 22 deg C:

`rate_scale = 1 / Q10^((37 - 22) / 10)`

This is a sensitivity analysis, not a claim that PyRhO's ChR2 fit was measured at 37 deg C.

## Calibration

For each Q10 condition:

1. Recalibrate `gbar` against Wang's reported saturated maximum current `Imax = 0.642 nA`.
2. Use saturating light `100 mW/mm2`, 100 ms, matching the previous calibration convention.
3. With that `gbar`, simulate Wang's 9.2 mW/mm2 1-s step.
4. Simulate Wang's seven intensity points: `0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2 mW/mm2`.

## Fixed recording/cell assumptions

- Cell: Dura-Bernal PT5B full morphology.
- Biophysics: original Dura-Bernal mechanisms.
- Slice: Wang-like in vitro slice abstraction.
- Optical attenuation: `mu_eff = 1.3 mm^-1`.
- Clamp: SEClamp at -70 mV.
- Rs: 10 MOhm.
- Electrode/amplifier readout: same post-filter as previous Wang suite.
- Temperature: 22 deg C.

