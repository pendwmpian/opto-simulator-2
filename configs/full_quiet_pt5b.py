"""Run-only overrides for the Dura-Bernal et al. quiet M1 network."""

QUIET_RATES_LONG = {
    "TPO": [0, 5], "TVL": [0, 2.5], "S1": [0, 5], "S2": [0, 5],
    "cM1": [0, 2.5], "M2": [0, 2.5], "OC": [0, 5],
}

FIXED_SEEDS = {"conn": 4321, "stim": 1234, "loc": 4321}
LFP_ELECTRODES_UM = [[150.0, 600.0, 150.0], [150.0, 800.0, 150.0], [150.0, 1000.0, 150.0]]


def apply_run_only_overrides(
    cfg, *, duration_ms, trial, output_dir, seeds=None, record_step_ms=0.1, record_lfp=True
):
    cfg.duration = float(duration_ms)
    cfg.seeds = dict(FIXED_SEEDS if seeds is None else seeds)
    cfg.ihGbar = 0.75
    cfg.makeKgbarFactorEqualToNewFactor = False
    cfg.ratesLong = {key: list(value) for key, value in QUIET_RATES_LONG.items()}
    cfg.addConn = 1
    cfg.addSubConn = 1
    cfg.addLongConn = 1
    cfg.addPulses = 0
    cfg.addIClamp = 0
    cfg.addNetStim = 0
    if hasattr(cfg, "pulse") and isinstance(cfg.pulse, dict):
        cfg.pulse["pop"] = "None"
    if hasattr(cfg, "pulse2") and isinstance(cfg.pulse2, dict):
        cfg.pulse2["pop"] = "None"
    cfg.saveCellSecs = False
    cfg.saveCellConns = False
    cfg.saveDataInclude = []
    cfg.savePickle = False
    cfg.saveJson = False
    cfg.gatherOnlySimData = True
    cfg.recordCells = []
    cfg.recordTraces = {
        "V_soma": {"sec": "soma", "loc": 0.5, "var": "v", "conds": {"pop": "PT5B"}}
    }
    cfg.recordStep = float(record_step_ms)
    cfg.recordLFP = [list(site) for site in LFP_ELECTRODES_UM] if record_lfp else []
    cfg.saveLFPPops = False
    cfg.saveLFPCells = False
    cfg.saveIMembrane = False
    cfg.recordDipoles = False
    cfg.recordDipole = False
    cfg.recordStim = False
    cfg.recordTime = False
    cfg.analysis = {}
    cfg.compactConnFormat = ["preGid", "sec", "loc", "synMech", "weight", "delay"]
    cfg.simLabel = f"quiet_pt5b_trial{trial}"
    cfg.saveFolder = str(output_dir)
    cfg.printRunTime = 0.1
    cfg.verbose = 0
    return cfg
