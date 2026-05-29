# Wang 2007 Fig. 3 Current-Clamp Initial Result

## Purpose

Test whether the calibrated Wang Fig. 2 photocurrent model can drive Fig. 3-like current-clamp spiking in the Dura-Bernal PT5B full cell.

## Model

```text
cell: PT5B_full
ChR2 gbar: 0.13478 mS/cm2
ChR2 kinetics: literature 4-state values
temperature: 22 degC
slice optical attenuation: mu_eff = 1.3 mm^-1
soma depth: 100 um
synaptic input: none
recording mode: current clamp
```

Voltage-clamp electrode/amplifier current filtering is not reused for spike detection. The biological readout is raw soma membrane voltage.

Script:

```text
scripts/run_wang_fig3_current_clamp.py
```

Outputs:

```text
outputs/wang2007_fig3_current_clamp/
outputs/wang2007_fig3_current_clamp_bias_0p2/
outputs/wang2007_fig3_current_clamp_bias_0p5/
outputs/wang2007_fig3_current_clamp_bias_1p0/
```

## Bug Found In Initial No-Spike Result

The first Fig. 3 run used `CellFromNetPyNE` in its older simplified mode. That loader did not apply the original Dura-Bernal/NetPyNE mechanisms stored in `PT5B_full_cellParams.pkl`. Instead it replaced the cell with a simplified `pas + hh` model and set leak reversal values near -70 mV.

This was the main reason the first current-clamp run did not spike. It was not caused by a leftover `SEClamp`.

The actual PT5B cellParams contain mechanisms such as:

```text
nax, kdr, kap, kdmc, hd, cal, can, kBK, cadad, pas, savedist
```

The loader now supports:

```text
use_original_biophysics=True
```

and the Fig. 3 script exposes:

```text
--use-original-biophysics
```

## Superseded Simplified-Biophysics Result

Protocol:

```text
1 s light flashes
intensity values:
  0.07, 0.14, 0.29, 0.58, 1.15, 2.3, 9.2 mW/mm2
DC bias current:
  0 nA
```

Result:

```text
resting Vm after equilibration:
  -72.1 mV

spike count:
  0 at all tested intensities

maximum Vm:
  0.07 mW/mm2 -> -68.6 mV
  0.14 mW/mm2 -> -67.1 mV
  0.29 mW/mm2 -> -65.4 mV
  0.58 mW/mm2 -> -63.9 mV
  1.15 mW/mm2 -> -62.6 mV
  2.3 mW/mm2 -> -61.0 mV
  9.2 mW/mm2 -> -57.2 mV
```

Interpretation:

The calibrated ChR2 current depolarizes the cell in an intensity-dependent way, but this PT5B model does not reach spike threshold from its natural resting state.

## DC Bias Sensitivity

Because Wang Fig. 3 is current clamp and baseline excitability is critical, small DC depolarizing bias was tested as a separate sensitivity analysis.

### 0.2 nA Bias

```text
rest Vm: -66.5 mV
spike count: 0 at all intensities
maximum Vm at 9.2 mW/mm2: -58.5 mV
```

### 0.5 nA Bias

```text
rest Vm: -62.6 mV
spike count: 0 at all intensities
maximum Vm at 9.2 mW/mm2: -59.0 mV
```

### 1.0 nA Bias

```text
rest Vm: -58.8 mV
spike count: 0 at all intensities
maximum Vm at 9.2 mW/mm2: -58.8 mV
```

Interpretation:

The lack of spiking is not solved by modest somatic depolarizing bias. This suggests that the issue is not only resting Vm, but also the current-clamp excitability/state of the imported PT5B model under this initialization and temperature.

## Frequency-Following Result

The first no-bias full run also attempted 4 ms pulse trains at:

```text
5, 10, 20, 30, 35, 40, 50 Hz
```

Result:

```text
evoked spike probability: 0 at all tested frequencies
```

This is expected given the no-bias 1 s flash result.

## Corrected Original-Biophysics Result

Protocol:

```text
cell: PT5B_full
biophysics: original Dura-Bernal/NetPyNE cellParams mechanisms
ChR2 gbar: 0.13478 mS/cm2
mu_eff: 1.3 mm^-1
temperature: 22 degC
DC bias: 0 nA
```

Output:

```text
outputs/wang2007_fig3_current_clamp_original_biophys/
outputs/wang2007_fig3_current_clamp_original_biophys_freq/
```

Intensity series:

```text
resting Vm after equilibration:
  -76.9 mV

0.07 mW/mm2:
  spike count = 0
  max Vm = -70.9 mV

0.14 mW/mm2:
  spike count = 0
  max Vm = -67.3 mV

0.29 mW/mm2:
  spike count = 0
  max Vm = -61.9 mV

0.58 mW/mm2:
  spike count = 1
  first AP peak latency = 26.8 ms
  max Vm = 37.5 mV

1.15 mW/mm2:
  spike count = 1
  first AP peak latency = 16.4 ms
  max Vm = 37.8 mV

2.3 mW/mm2:
  spike count = 1
  first AP peak latency = 11.8 ms
  max Vm = 37.8 mV

9.2 mW/mm2:
  spike count = 1
  first AP peak latency = 7.8 ms
  max Vm = 37.9 mV
```

Hill fit to spike count is poorly constrained because the count saturates at one spike:

```text
K ~= 0.42 mW/mm2
Imax ~= 1.03 spikes
Hill n hits upper fit bound
```

Frequency following with 4 ms pulses at 9.2 mW/mm2:

```text
5 Hz:
  success probability = 0.40

10 Hz:
  success probability = 0.20

20 Hz:
  success probability = 0.10

30 Hz:
  success probability = 0.067

35 Hz:
  success probability = 0.057

40 Hz:
  success probability = 0.05

50 Hz:
  success probability = 0.04
```

## Corrected Interpretation

With original Dura-Bernal biophysics, the Fig. 2-calibrated ChR2 conductance can drive spikes in the PT5B model. The previous no-spike result should be discarded.

The corrected model captures several qualitative Wang Fig. 3 features:

```text
spike threshold occurs around sub-mW/mm2 irradiance
first AP latency decreases with light intensity
9.2 mW/mm2 first AP peak latency is several ms
```

Current mismatches:

```text
Wang maximum AP count during 1 s flashes:
  about 25 spikes

current PT5B model:
  one spike only across suprathreshold intensities

Wang frequency following:
  reliable to about 30 Hz

current PT5B model:
  poor following; only about 2 evoked spikes per train
```

This points to intrinsic current-clamp excitability/adaptation of the PT5B model, not ChR2 current scale.

## Updated Recommended Next Step

The next diagnostic remains direct current injection, but now for a different reason:

```text
measure rheobase
measure f-I curve
check whether the original PT5B model can fire repetitively at 22 degC
```

If direct current injection also produces only one spike, then Fig. 3 mismatch is dominated by intrinsic adaptation/excitability of the chosen PT5B model.

If direct current injection produces repetitive firing but ChR2 does not, then the issue is the spatial distribution and time course of optogenetic current.

## Superseded Interpretation

The Fig. 2-calibrated ChR2 conductance is not sufficient to reproduce Wang Fig. 3 spiking in the current PT5B model as initialized here.

This should not yet be interpreted as a failure of the ChR2 `gbar` calibration. Fig. 2 voltage-clamp photocurrent calibration is working well. The problem is likely current-clamp excitability.

Likely contributors:

```text
1. Dura-Bernal PT5B model excitability differs from Wang's recorded L5 pyramidal neurons.
2. Initial Na/K/HCN channel states may not be appropriate for current-clamp spike tests at 22 degC.
3. The cell may require rheobase/current-injection calibration before optogenetic spiking validation.
4. Active channel conductances from the network model may not be tuned for acute-slice current clamp.
5. Somatic DC bias alone may be insufficient if dendritic/axonal spike initiation is not represented appropriately.
```

## Recommended Next Step

Before using Fig. 3 as an optogenetic validation, run a current-injection excitability check:

```text
inject somatic square current pulses:
  0.1-2.0 nA

measure:
  rheobase
  spike threshold
  first-spike latency
  firing frequency
  whether the model spikes at all at 22 degC
```

If the cell cannot spike robustly to direct current injection, Fig. 3 cannot be meaningfully simulated without first fixing or selecting a more appropriate current-clamp cell model.

If the cell can spike to direct current but not to ChR2, then the problem is optogenetic current distribution/amplitude rather than intrinsic excitability.
