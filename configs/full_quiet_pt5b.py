"""Run-only overrides for the Dura-Bernal et al. quiet M1 network."""

QUIET_RATES_LONG = {
    "TPO": [0, 5], "TVL": [0, 2.5], "S1": [0, 5], "S2": [0, 5],
    "cM1": [0, 2.5], "M2": [0, 2.5], "OC": [0, 5],
}


def apply_run_only_overrides(cfg, *, duration_ms, trial, output_dir):
    seed_offset = 17 * trial
    cfg.duration = float(duration_ms)
    cfg.seeds = {"conn": 4321 + seed_offset, "stim": 1234 + seed_offset, "loc": 4321 + seed_offset}
    cfg.ratesLong = {key: list(value) for key, value in QUIET_RATES_LONG.items()}
    cfg.addPulses = 0
    cfg.addIClamp = 0
    cfg.addNetStim = 0
    cfg.saveCellSecs = False
    cfg.saveCellConns = False
    cfg.saveDataInclude = ["simData", "simConfig"]
    cfg.savePickle = False
    cfg.saveJson = False
    cfg.gatherOnlySimData = True
    cfg.recordCells = []
    cfg.recordTraces = {"V_soma": {"sec": "soma", "loc": 0.5, "var": "v"}}
    cfg.recordStep = 0.1
    cfg.recordLFP = []
    cfg.saveLFPPops = False
    cfg.saveLFPCells = False
    cfg.recordDipoles = False
    cfg.recordStim = False
    cfg.recordTime = False
    cfg.analysis = {}
    cfg.compactConnFormat = ["preGid", "sec", "loc", "synMech", "weight", "delay"]
    cfg.simLabel = f"quiet_pt5b_trial{trial}"
    cfg.saveFolder = str(output_dir)
    cfg.printRunTime = 0.1
    cfg.verbose = 0
    return cfg
