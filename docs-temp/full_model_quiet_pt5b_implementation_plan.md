# Full Model Quiet PT5B Baseline Implementation Plan

## 1. Goal

Dura-Bernal et al. (2023) の NetPyNE/NEURON full M1 network を維持し、まず通常・quiet 状態で 1 秒の simulation を最大 3 trial 実行する。主な観測対象は以下に限定する。

- 全細胞の spike time / gid
- あらかじめ固定した PT5B target cell 1 個の soma Vm
- wall time、CPU time、peak RSS、rank 数、終了状態

第一段階では ChR2 を発現・照射しない baseline を作り、Schiemann et al. (2015) と Dura-Bernal et al. (2023) が用いた quiet-state PT5B firing の再現性を評価する。第二段階で、同じ network baseline に Thy1-ChR2 発現と光刺激を加える。

## 2. Current Repository Audit

2026-06-12 時点のワークツリーには `docs/`、`scripts/`、`mod/` はあるが、既存文書が参照する次の実体が存在しない。

- `external/M1_NetPyNE_CellReports_2023/`
- `pyproject.toml`
- `uv.lock`
- full model 用に compile 済みの NMODL mechanisms

したがって、実装は upstream model の再取得、commit 固定、Python/NEURON/NetPyNE/MPI 環境の固定から始める。既存の単一細胞スクリプトが参照する path は、そのままでは現在実行できない。

Upstream README は full model について 96 cores で 1 秒あたり約 2 時間と記載している。8 cores では単純換算でも約 24 時間だが、通信・memory bandwidth・load imbalance によりさらに長くなる可能性がある。最初から 3 trial を投入せず、smoke test と 1 trial の実測後に判断する。

## 3. Scientific Baseline

quiet baseline は upstream の Figure 2 "Control Quiet" 設定を出発点にする。upstream `cfg.py` の主な既定値は次の通り。

```python
cfg.duration = 1000.0
cfg.dt = 0.025
cfg.hParams = {'celsius': 34, 'v_init': -80}
cfg.scale = 0.3
cfg.cellmod['PT5B'] = 'HH_full'
cfg.addConn = 1
cfg.addSubConn = 1
cfg.addLongConn = 1
cfg.addPulses = 1
cfg.ratesLong = {
    'TPO': [0, 5], 'TVL': [0, 2.5], 'S1': [0, 5],
    'S2': [0, 5], 'cM1': [0, 2.5], 'M2': [0, 2.5], 'OC': [0, 5],
}
```

実装時には Figure 2 quiet batch の全 override を upstream commit から抽出し、値と出典行を manifest に保存する。baseline では movement pulse と ChR2 irradiance を 0 にし、ChR2 mechanism が network に組み込まれていても暗状態で電流を流さないことを確認する。

### Trial definition

- trial 0: upstream quiet seed をそのまま使用する再現 trial
- trial 1-2: `conn`、`stim`、`loc` seed を明示的に変えた biological variability trial
- target gid: trial 間で同一 population 内の同じ deterministic selection rule により決定する
- target selection: PT5B の local gid 一覧から中央に近い cell を 1 個選ぶことを第一候補とし、選択規則、gid、position を manifest に記録する

## 4. Run-Only Minimal Configuration Review

提案された削減方針は採用する。ただし NetPyNE version による属性名の差を実行前 test で確認する。特に upstream は `recordDipoles`、提案は `recordDipole` であり、実際にインストールした version の API を正とする。

```python
cfg.saveCellSecs = False
cfg.saveCellConns = False
cfg.saveDataInclude = ['simData', 'simConfig']

cfg.recordCells = [target_gid]
cfg.recordTraces = {
    'V_soma': {'sec': 'soma', 'loc': 0.5, 'var': 'v'},
}
cfg.recordStep = 0.1

cfg.recordLFP = []
cfg.saveLFPPops = False
cfg.saveLFPCells = False
cfg.recordDipoles = False
cfg.recordStim = False
cfg.recordTime = False
cfg.analysis = {}

cfg.compactConnFormat = [
    'preGid', 'sec', 'loc', 'synMech', 'weight', 'delay'
]
```

`recordStep=0.1 ms` では 1 秒あたり約 10,001 点だけなので、target 1 cell に限定すれば小さい。Vm waveform の fidelity に問題がなければ 0.2 ms を第二候補とする。1.0 ms は spike waveform と threshold crossing の評価には粗いため、OOM 対策ではなく最終 fallback とする。

### Important memory distinction

- `saveCellSecs=False`、`saveCellConns=False`、`saveDataInclude` は主に gather/save 時の出力を減らす。
- `recordCells`、LFP/dipole/analysis 無効化は recording buffer と解析時 memory を減らす。
- `compactConnFormat` は connection の Python metadata を小さくできる可能性があり、`connectCells()` 中の memory に直接関係する。
- NetCon、synapse mechanism、cell sections は simulation に必要なので削除しない。
- `createPyStruct=False` は大きな削減候補だが、NetPyNE の connection construction、subcellular placement、recording、saving との互換性を version ごとに確認するまで採用しない。

## 5. OOM Avoidance Strategy

OOM は save 時ではなく `sim.net.connectCells()` 中に起こり得るため、network construction を段階計測する専用 runner を作る。

1. `initialize`
2. `createPops`
3. `createCells`
4. `connectCells`
5. `addStims`
6. `setupRecording`
7. `runSim`
8. spike/Vm の minimal gather と save

各段階の前後で rank ごとの RSS と host 全体の available memory を記録する。rank 0 の RSS だけでなく、全 rank の合計と最大値を report する。

採用順序は以下とする。

1. 8 MPI ranks で cell/connection ownership を分散する。
2. `compactConnFormat` を connection creation 前から有効にする。
3. `saveCellSecs=False`、`saveCellConns=False`、minimal recording を有効にする。
4. built-in distributed node saving が不要な metadata を作る場合は、spike と target trace だけを明示的に gather/save する runner に置き換える。
5. `createPyStruct=False` または connection Python metadata の post-build release は、10-100 ms small-scale equivalence test で spike/Vm が一致した場合のみ採用する。

full model、cell density、morphology、ion channels、synMech、NetCon、background input は baseline の科学的条件なので、OOM 回避のために scale や synapse 数を黙って減らさない。scale reduction は pipeline smoke test と明記した場合だけ使用する。

## 6. Implementation Artifacts

以下を追加する予定。

```text
external/M1_NetPyNE_CellReports_2023/   pinned upstream checkout
pyproject.toml                           uv-managed environment
uv.lock                                  exact Python dependencies
scripts/run_full_quiet_pt5b.py           staged NetPyNE runner
scripts/monitor_full_run.py               periodic process/memory sampler
scripts/summarize_quiet_pt5b.py           minimal scientific summary
configs/full_quiet_pt5b.py                run-only SimConfig overrides
docs-temp/full-model-runs/<run_id>/
```

`run_id` は condition、trial、seed、timestamp を含む。各 run directory に以下だけを置く。

```text
manifest.json
sim_config.json
spikes.npz
target_vm.npz
runtime.json
resource_usage.csv
stdout.log
stderr.log
summary.md
```

full `net`、`netParams`、全 cell section、全 connection metadata、LFP、dipole、stimulus trace、plot object は保存しない。再現に必要な upstream commit、local patch hash、software versions、config overrides、seed、target gid、command line は `manifest.json` に保存する。

## 7. Environment and Compatibility Work

1. `uv venv --python=3.13` を作成する。
2. `uv add` で NetPyNE、NEURON、NumPy、SciPy と必要最小限の解析依存を追加する。
3. Python 3.13 と upstream-era NetPyNE の互換性を確認する。非互換なら、Python version を勝手に変更せず、理由と compatible version を記録して環境方針を更新する。
4. MPI-enabled NEURON と `mpiexec` の smoke test を行う。
5. upstream NMODL と local ChR2 MOD を同じ mechanism build に含め、duplicate mechanism 名を検査する。
6. NetPyNE version ごとに `compactConnFormat`、dipole/LFP flags、recording selection、minimal save の挙動を unit/smoke test する。

## 8. Execution and pueue Monitoring

long run は `pueue` に投入する。実装後の command 形は概ね以下とする。

```bash
rtk proxy pueue add 'systemd-run --user --scope -p MemoryMax=220G -p MemoryHigh=200G /usr/bin/time -v mpiexec -np 8 uv run python scripts/run_full_quiet_pt5b.py --trial 0 --duration-ms 1000'
```

実環境で user scope が利用できない場合は、`ulimit`/MPI launcher の利用可能な resource control を調査し、同等の hard cap を runner wrapper に実装する。初期上限は 220 GB、soft pressure threshold は 200 GB とし、OS と Codex server 用に約 36 GB を残す。CPU は `mpiexec -np 8` で 8 ranks に固定し、BLAS/OpenMP の oversubscription を防ぐため `OMP_NUM_THREADS=1` 等も固定する。

監視項目:

- pueue task id / state / exit code
- wall clock、user/system CPU、maximum resident set size
- 30-60 秒間隔の total RSS、rank 最大 RSS、system available memory
- simulation phase と NEURON simulation time
- output directory size
- timeout、MPI abort、OOM kill の検出

`pueue status` と `pueue log` は run manifest に対応付ける。監視 script は memory threshold 超過時に新規 trial を開始せず、実行中 job は hard cap によって Codex server を保護する。

## 9. Staged Run Plan

### Stage A: Static and API validation

- upstream commit と Figure 2 quiet config を固定する
- config attribute compatibility を確認する
- target PT5B selection を確認する
- ChR2 mechanism が dark baseline で zero current になることを確認する

### Stage B: Construction-only smoke tests

- reduced scale、10 ms、1-2 ranks で phase instrumentation を検証する
- full scale、simulation なしで `connectCells()` 完了まで 8 ranks で実行する
- peak memory が 200 GB を超える見込みなら trial 本番へ進まない

### Stage C: Short full-network dynamics

- full scale、50-100 ms、8 ranks
- 全 spike と target soma Vm が正しく保存されることを確認する
- minimal output と通常 gather の small test が一致することを確認する

### Stage D: One-second quiet baseline

- trial 0 を 1 秒実行する
- runtime と peak memory を実測する
- scientific acceptance と resource budget を満たした場合だけ trial 1-2 を逐次投入する
- 3 trial を同時実行しない

### Stage E: ChR2 extension

- quiet baseline の network、seed handling、recording を固定したまま ChR2 expression と light stimulation だけを追加する
- baseline dark-current equivalence test 後に光条件へ進む

## 10. Evaluation Criteria

### Technical acceptance

- `connectCells()` が OOM なしで完了する
- 8 ranks の simulation が正常終了する
- output に全細胞 spike と target PT5B soma Vm が含まれる
- full network metadata と不要 trace が output に含まれない
- run manifest から commit、config、seed、command、target gid を再構成できる
- runtime と peak memory が記録される

### Quiet PT5B scientific acceptance

各 trial と trial aggregate について以下を算出する。

- PT5B population firing-rate distribution
- mean、median、SD、IQR、silent-cell fraction
- target PT5B firing rate、ISI、burst count
- target soma Vm の mean、SD、subthreshold range、spike count
- 全 population の mean firing rate
- transient を除外した評価 window と、1 秒全体の両方

Schiemann et al. (2015) の experimental distribution と Dura-Bernal et al. (2023) Figure 2 の digitized/quoted target を、実装前に一次資料から表形式で固定する。1 秒・最大 3 trial では low-rate PT5B の統計誤差が大きいため、厳密な同等性ではなく以下を判定する。

- population mean/median が experimental range と同じ order にあるか
- silent fraction と firing-rate distribution の形が大きく矛盾しないか
- target cell が quiet network 内で非生理的な持続 depolarization、異常 burst、numerical instability を示さないか
- trial 間差が target mismatch より大きくないか

Schiemann target 値を未確認のまま acceptance threshold を推測して固定しない。

## 11. Reporting

最終的な `docs-temp` package には以下をまとめる。

- 実行目的と再現対象
- upstream/local commit と software versions
- full config override と baseline assumptions
- target cell selection
- trial ごとの seed、runtime、peak memory、exit status
- PT5B および全 population spike statistics
- target Vm summary と必要最小限の figure
- Schiemann/Dura-Bernal target との比較
- deviations、known limitations、OOM mitigation の効果
- ChR2 stage へ引き継ぐ固定条件

figure は simulation 中には作らず、保存済みの `spikes.npz` と `target_vm.npz` から post hoc に生成する。`docs-temp/full-model-runs/` 以下だけを archive すれば、結果と仮定をまとめてダウンロードできる構成にする。

## 12. Main Risks and Decisions

- **Runtime risk:** upstream estimate から、8 cores で 1 秒は 1 日以上かかる可能性がある。trial 数は実測後に 1-3 で決める。
- **Connection-memory risk:** save flags だけでは construction OOM を解決しない。MPI distribution と compact connection representation を先に検証する。
- **Version risk:** upstream model は古い NetPyNE API を前提としている可能性がある。最新 version へ無条件移植せず、動作する version を lock する。
- **Statistical risk:** 1 秒 x 3 trial は quiet low-rate cells の精密推定には短い。初期 feasibility baseline として扱う。
- **Scientific scope:** full network の cell model、morphology、ion channels、synapses、NetCon、background input は保持する。network scale や synapse count の削減は本番結果には採用しない。

## 13. References Used for This Plan

- Dura-Bernal et al. model repository: <https://github.com/suny-downstate-medical-center/M1_NetPyNE_CellReports_2023>
- Upstream `cfg.py`: <https://github.com/suny-downstate-medical-center/M1_NetPyNE_CellReports_2023/blob/main/sim/cfg.py>
- Upstream `batch.py`: <https://github.com/suny-downstate-medical-center/M1_NetPyNE_CellReports_2023/blob/main/sim/batch.py>
- Local prior work: `docs/phase0_dura_bernal_model.md`, `docs/walkthrough.md`, and the existing ChR2/Wang validation documents

Schiemann et al. (2015) の exact bibliographic entry、quiet PT5B sample definition、numerical target は実装 Stage A で一次資料から確定し、report に DOI/PMID とともに追記する。
