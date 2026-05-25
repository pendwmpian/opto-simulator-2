# Conductance-Based ChR2 Model

## Concept

The current exploratory model injects a fixed inward current:

```text
I = constant(local irradiance, expression)
```

This is useful for spatial pipeline testing, but it is not a real ion channel model.

A conductance-based ChR2 model instead uses:

```text
I_ChR2 = g_ChR2(t, light, expression) * (V - E_ChR2)
```

or equivalently, with inward-current sign convention:

```text
I_inward = g_ChR2(t, light, expression) * (E_ChR2 - V)
```

This matters because the current depends on membrane voltage. When the membrane depolarizes toward `E_ChR2`, the driving force becomes smaller.

## What Is `gbar`?

`gbar` means maximum conductance.

Conductance is the inverse of resistance:

```text
conductance = 1 / resistance
```

In the ChR2 context, `gbar` summarizes how much current can pass through all ChR2 channels in a membrane patch when the channels are maximally available/open.

Biologically, `gbar` absorbs:

- number of ChR2 molecules
- single-channel conductance
- membrane expression density
- trafficking/localization
- scaling mismatch between model and real cells

For this project, `gbar` is the main fitted parameter for ChR2 expression/sensitivity.

## Why Is `E_ChR2` Around 0 mV?

`E_ChR2` is the reversal potential of the ChR2 current.

At this membrane voltage, net ChR2 current is zero:

```text
I_ChR2 = g * (V - E_ChR2) = 0 when V = E_ChR2
```

ChR2 is a nonselective cation channel. It passes mainly positive ions such as Na+, K+, H+, and Ca2+. Because it is not a pure Na or pure K channel, its reversal potential lies between the major cation equilibrium potentials and is commonly approximated near 0 mV.

Interpretation:

```text
V = -70 mV: E_ChR2 - V ~= 70 mV, strong inward depolarizing current
V = -20 mV: E_ChR2 - V ~= 20 mV, weaker inward current
V = 0 mV: no net ChR2 current
V > 0 mV: current can reverse outward
```

This is why fixed-current injection overestimates ChR2 current during strong depolarization and spikes.

## Literature Basis

### Nikolic et al., 2009

Nikolic et al. introduced a ChR2 photocycle model with closed and open states, including adaptation/desensitization dynamics. The key point for this project is that ChR2 current cannot be represented only by instantaneous light intensity; channel state history matters.

Reference:

```text
Nikolic K, Grossman N, Grubb MS, Burrone J, Toumazou C, Degenaar P.
Photocycles of channelrhodopsin-2.
Photochemistry and Photobiology, 2009.
```

### Foutz et al., 2012

Foutz et al. coupled a ChR2 membrane dynamics model to a multicompartment cortical pyramidal neuron and an optical model. Their result is directly relevant: stimulation threshold depends strongly on illuminated membrane surface area and source proximity.

They model ChR2 as a nonspecific cation channel with reversal potential near 0 mV and use channel states to compute photocurrent.

Reference:

```text
Foutz TJ, Arlow RL, McIntyre CC.
Theoretical principles underlying optical stimulation of a channelrhodopsin-2 positive pyramidal neuron.
Journal of Neurophysiology, 2012.
```

### Grossman et al., 2013

Grossman et al. implemented ChR2 in NEURON and studied how spatial illumination pattern changes spike kinetics. This is very close to the current project.

They emphasize that maximum conductance accounts for both single-channel conductance and expression level, which matches our use of `gbar` as the fit parameter.

Reference:

```text
Grossman N, Simiaki V, Martinet C, Toumazou C, Schultz SR, Nikolic K.
The spatial pattern of light determines the kinetics and modulates backpropagation of optogenetic action potentials.
Journal of Computational Neuroscience, 2013.
```

### Foutz et al., 2014

Foutz et al. give a clear four-state ChR2 implementation for optogenetic stimulation models, including `E_cat = 0 mV` and O1/O2 channel conductances.

Reference:

```text
Foutz TJ, Arlow RL, McIntyre CC.
Theoretical principles underlying optical stimulation of myelinated axons expressing channelrhodopsin-2.
```

## Implementation Choice

The next implementation should use a compact four-state conductance model:

```text
C1 <-> O1 <-> O2 <-> C2 <-> C1
```

For the first implementation:

- optics are fixed from Yona/near-488 nm assumptions
- kinetics are fixed from the ChR2 literature
- `gbar` is the fitted ChR2 expression/sensitivity parameter
- `E_ChR2 = 0 mV`

This is a good balance between biological realism and computational cost.

