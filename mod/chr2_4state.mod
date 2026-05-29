NEURON {
    SUFFIX chr2_4state
    NONSPECIFIC_CURRENT i
    RANGE gbar, e, gamma, irr
    RANGE k1, k2, gd1, gd2, e12, e21, e12_dark, e21_dark, gr
    RANGE q10_scale
    RANGE q10_k1_scale, q10_k2_scale, q10_gd1_scale, q10_gd2_scale
    RANGE q10_e12_scale, q10_e21_scale, q10_gr_scale
    RANGE fit_e12_scale, fit_e21_scale, fit_gd2_scale, fit_gamma_scale
    RANGE activation_mode, wavelength_nm, wloss, sigma_retinal, tau_chr2
    RANGE epsilon1, epsilon2
    RANGE c1, o1, o2, c2, p, open
    RANGE c1_init, o1_init, o2_init, c2_init
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mW) = (milliwatt)
}

PARAMETER {
    gbar = 0.0 (millimho/cm2)
    e = 0.0 (mV)
    gamma = 0.05
    irr = 0.0 (mW/mm2)

    k1 = 0.35 (/ms) : scaffold only, ignored in photon activation mode
    k2 = 0.12 (/ms) : scaffold only, ignored in photon activation mode
    gd1 = 0.130 (/ms)
    gd2 = 0.025 (/ms)
    e12 = 0.053 (/ms)
    e21 = 0.023 (/ms)
    e12_dark = 0.022 (/ms)
    e21_dark = 0.011 (/ms)
    gr = 0.0004 (/ms)
    activation_mode = 0
    wavelength_nm = 470
    wloss = 1.3
    sigma_retinal = 12e-20 (m2)
    tau_chr2 = 1.3 (ms)
    epsilon1 = 0.8535
    epsilon2 = 0.14
    c1_init = 1.0
    o1_init = 0.0
    o2_init = 0.0
    c2_init = 0.0
    q10_scale = 1.0
    q10_k1_scale = 1.0
    q10_k2_scale = 1.0
    q10_gd1_scale = 1.0
    q10_gd2_scale = 1.0
    q10_e12_scale = 1.0
    q10_e21_scale = 1.0
    q10_gr_scale = 1.0
    fit_e12_scale = 1.0
    fit_e21_scale = 1.0
    fit_gd2_scale = 1.0
    fit_gamma_scale = 1.0
}

ASSIGNED {
    v (mV)
    i (mA/cm2)
    open
    ka1 (/ms)
    ka2 (/ms)
    e12_eff (/ms)
    e21_eff (/ms)
    ephoton (joule)
    flux (/m2-s)
    fphi (/ms)
    s0
    hc (joule-m)
}

STATE {
    c1
    o1
    o2
    c2
    p
}

INITIAL {
    c1 = c1_init
    o1 = o1_init
    o2 = o2_init
    c2 = c2_init
    p = 0.0
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    open = o1 + fit_gamma_scale * gamma * o2
    i = (0.001 * gbar) * open * (v - e)
}

DERIVATIVE states {
    ka1 = q10_scale * q10_k1_scale * k1 * irr
    ka2 = q10_scale * q10_k2_scale * k2 * irr
    s0 = 0.0
    if (activation_mode > 0.5) {
        hc = 1.986446e-25
        ephoton = 1e9 * hc / wavelength_nm
        flux = 1000.0 * irr / ephoton
        fphi = flux * sigma_retinal / (wloss * 1000.0)
        s0 = 0.5 * (1.0 + tanh(120.0 * (100.0 * irr - 0.1)))
        ka1 = q10_scale * q10_k1_scale * epsilon1 * fphi * p
        ka2 = q10_scale * q10_k2_scale * epsilon2 * fphi * p
    }
    p' = (s0 - p) / tau_chr2

    e12_eff = e12
    e21_eff = e21
    if (irr <= 0.0) {
        e12_eff = e12_dark
        e21_eff = e21_dark
    }

    c1' = -ka1 * c1 + q10_scale * q10_gd1_scale * gd1 * o1 + q10_scale * q10_gr_scale * gr * c2
    o1' = ka1 * c1 - q10_scale * q10_gd1_scale * gd1 * o1 - q10_scale * q10_e12_scale * fit_e12_scale * e12_eff * o1 + q10_scale * q10_e21_scale * fit_e21_scale * e21_eff * o2
    o2' = ka2 * c2 + q10_scale * q10_e12_scale * fit_e12_scale * e12_eff * o1 - q10_scale * q10_e21_scale * fit_e21_scale * e21_eff * o2 - q10_scale * q10_gd2_scale * fit_gd2_scale * gd2 * o2
    c2' = q10_scale * q10_gd2_scale * fit_gd2_scale * gd2 * o2 - ka2 * c2 - q10_scale * q10_gr_scale * gr * c2
}
