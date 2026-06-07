NEURON {
    SUFFIX chr2_pyrho6
    NONSPECIFIC_CURRENT i
    RANGE gbar, e, gamma, irr, wavelength_nm
    RANGE phi_m, phi_m_scale, k1, k2, p_exp, gf0, kf, gb0, kb, q_exp, go1, go2, gd1, gd2, gr0
    RANGE q10_scale
    RANGE c1, i1, o1, o2, i2, c2, open, flux
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mW) = (milliwatt)
}

PARAMETER {
    gbar = 0.0 (millimho/cm2)
    e = 0.0 (mV)
    gamma = 0.00369
    irr = 0.0 (mW/mm2)
    wavelength_nm = 470

    : PyRhO six-state ChR2 fit values, units are ms^-1 except phi_m.
    : Table values follow Evans et al. 2016 / PyRhO ChR2 six-state example fit.
    phi_m = 5.02e17
    phi_m_scale = 1.0
    k1 = 18.2 (/ms)
    k2 = 4.07 (/ms)
    p_exp = 0.981
    gf0 = 0.0365 (/ms)
    kf = 0.121 (/ms)
    gb0 = 0.0143 (/ms)
    kb = 0.131 (/ms)
    q_exp = 1.45
    go1 = 1.93 (/ms)
    go2 = 3.38 (/ms)
    gd1 = 0.108 (/ms)
    gd2 = 0.0115 (/ms)
    gr0 = 0.00033 (/ms)

    q10_scale = 1.0
}

ASSIGNED {
    v (mV)
    i (mA/cm2)
    open
    flux
    ga1 (/ms)
    ga2 (/ms)
    gf (/ms)
    gb (/ms)
    hp
    hq
    ephoton (joule)
    hc (joule-m)
}

STATE {
    c1
    i1
    o1
    o2
    i2
    c2
}

INITIAL {
    c1 = 1.0
    i1 = 0.0
    o1 = 0.0
    o2 = 0.0
    i2 = 0.0
    c2 = 0.0
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    open = o1 + gamma * o2
    i = (0.001 * gbar) * open * (v - e)
}

DERIVATIVE states {
    hc = 1.986446e-25
    ephoton = 1e9 * hc / wavelength_nm
    flux = 1000.0 * irr / ephoton / 1.0e6 : photons/mm2/s

    hp = 0.0
    hq = 0.0
    if (flux > 0.0) {
        hp = flux^p_exp / (flux^p_exp + (phi_m * phi_m_scale)^p_exp)
        hq = flux^q_exp / (flux^q_exp + (phi_m * phi_m_scale)^q_exp)
    }

    ga1 = q10_scale * k1 * hp
    ga2 = q10_scale * k2 * hp
    gf = q10_scale * (kf * hq + gf0)
    gb = q10_scale * (kb * hq + gb0)

    c1' = q10_scale * gd1 * o1 + q10_scale * gr0 * c2 - ga1 * c1
    i1' = ga1 * c1 - q10_scale * go1 * i1
    o1' = q10_scale * go1 * i1 + gb * o2 - (q10_scale * gd1 + gf) * o1
    o2' = q10_scale * go2 * i2 + gf * o1 - (q10_scale * gd2 + gb) * o2
    i2' = ga2 * c2 - q10_scale * go2 * i2
    c2' = q10_scale * gd2 * o2 - (q10_scale * gr0 + ga2) * c2
}
