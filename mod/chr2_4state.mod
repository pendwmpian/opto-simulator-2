NEURON {
    SUFFIX chr2_4state
    NONSPECIFIC_CURRENT i
    RANGE gbar, e, gamma, irr
    RANGE k1, k2, gd1, gd2, e12, e21, gr
    RANGE c1, o1, o2, c2, open
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mW) = (milliwatt)
}

PARAMETER {
    gbar = 0.0 (millimho/cm2)
    e = 0.0 (mV)
    gamma = 0.1
    irr = 0.0 (mW/mm2)

    k1 = 0.35 (/ms) : per mW/mm2
    k2 = 0.12 (/ms) : per mW/mm2
    gd1 = 0.10 (/ms)
    gd2 = 0.025 (/ms)
    e12 = 0.011 (/ms)
    e21 = 0.008 (/ms)
    gr = 0.00033 (/ms)
}

ASSIGNED {
    v (mV)
    i (mA/cm2)
    open
    ka1 (/ms)
    ka2 (/ms)
}

STATE {
    c1
    o1
    o2
    c2
}

INITIAL {
    c1 = 1.0
    o1 = 0.0
    o2 = 0.0
    c2 = 0.0
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    open = o1 + gamma * o2
    i = (0.001 * gbar) * open * (v - e)
}

DERIVATIVE states {
    ka1 = k1 * irr
    ka2 = k2 * irr

    c1' = -ka1 * c1 + gd1 * o1 + gr * c2
    o1' = ka1 * c1 - gd1 * o1 - e12 * o1 + e21 * o2
    o2' = ka2 * c2 + e12 * o1 - e21 * o2 - gd2 * o2
    c2' = gd2 * o2 - ka2 * c2 - gr * c2
}
