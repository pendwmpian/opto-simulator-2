# NEURON / CoreNEURON 詳細比較レポート

作成日: 2026-07-13  
対象ブランチ: `codex/coreneuron-backend-comparison`

## 1. 結論

CoreNEURON は、外部入力が弱く物理細胞が発火しない 100 ms 条件では、通常
NEURON と全スパイク・全 PT5B soma Vm がビット単位で一致した。この条件では
積分部分が 3.05 倍、総時間が 2.48 倍高速だった。

一方、3 s 実行と同じ入力系列を用いた発火ありの先頭 100 ms 条件では、入力、
seed、接続生成、MOD ソース、MPI 数を一致させても CoreNEURON のみが最初の
ネットワーク発火段階から分岐した。

- 外部入力: 両方とも 416 発、GID・時刻がビット単位で一致
- 通常 NEURON の物理細胞スパイク: 586 発
- CoreNEURON の物理細胞スパイク: 16,270 発
- 全 PT5B Vm の最初の差: 47.0 ms
- 100 ms 共通区間の全 PT5B Vm RMSE: 17.869 mV
- 最大絶対差: 100.628 mV

互換化済み MOD を通常 NEURON で実行した結果は、June 17 JST の NEURON
基準とスパイク・Vm ともビット単位で一致した。したがって、観測差は seed、
LC/NA 相当設定、MOD 互換化そのものではなく、モデルを CoreNEURON へ転送・
実行する経路に限定される。

現状の CoreNEURON 結果は、この発火ありモデルの科学的解析には使用できない。

## 2. 比較対象

### June NEURON 基準

- 出力: `docs-temp/full-model-runs/full-quiet-3s-ih075-16r-pt5ball-trial0`
- NEURON 9.0.1
- NetPyNE 1.1.1
- 16 MPI ranks
- 3,000 ms、`dt=0.025 ms`、記録間隔 `0.1 ms`
- LFP なし
- 全 PT5B 1,435 細胞の soma Vm を記録

### CoreNEURON 3 s 実行

- 出力: `docs-temp/full-model-runs/full-quiet-3s-coreneuron827-16r-trial0`
- Python パッケージ上の NEURON 8.2.7
- CoreNEURON 実行バナー: `8.2.2 dc3ead2f`
- NetPyNE 1.1.1
- CPU、16 MPI ranks、OpenMP 1 thread/rank
- GPU 無効、`cell-permute=0`
- 3,000 ms、`dt=0.025 ms`、記録間隔 `0.1 ms`
- LFP なし

### matched-prefix 対照

3 s 条件と同じ VecStim 系列を生成後、積分時間だけを 100 ms に戻して比較した。

- NEURON:
  `docs-temp/full-model-runs/full-quiet-prefix100ms-input3s-neuron827-patched-16r-trial0`
- 比較対象 CoreNEURON:
  3 s CoreNEURON 出力の先頭 100 ms

比較用ランナーは `--input-duration-ms` で入力生成時間幅と積分時間幅を分離する。

## 3. メモリ使用量

### 3.1 共通 phase の process RSS 合計

各 `resource_rank*.jsonl` に保存された process RSS を、同じ phase の 16 rank で
合計した。GiB は `2^30 bytes` で換算している。

| Phase | NEURON 合計 RSS | CoreNEURON 合計 RSS | CoreNEURON 平均/rank |
|---|---:|---:|---:|
| cell 作成後 | 12.42 GiB | 9.84 GiB | 0.62 GiB |
| 接続作成後 | 77.32 GiB | 67.34 GiB | 4.21 GiB |
| 記録設定後 | 77.63 GiB | 68.62 GiB | 4.29 GiB |
| simulation 後 | 80.61 GiB | 74.28 GiB | 4.64 GiB |
| gather 後 | 59.72 GiB | 59.46 GiB | 3.72 GiB |

共通の phase 境界で観測された CoreNEURON の最大同時 RSS は **74.28 GiB**。
通常 NEURON の 80.61 GiB より 6.33 GiB、7.86% 少なかった。記録設定直後では
9.01 GiB、11.61% 少なかった。

simulation 後の CoreNEURON は、rank あたり平均 4.64 GiB、最大 4.70 GiB だった。
gather 中は rank 0 が 8.93 GiB まで増えたが、他 rank が減少していたため、同じ
phase の合計は 59.46 GiB だった。

### 3.2 CoreNEURON 内部 telemetry の一時ピーク

CoreNEURON の内部ログは `After nrn_setup` で以下を報告した。

- rank 平均: 5,235.8933 MB
- rank 最小: 5,165.0352 MB
- rank 最大: 5,287.4805 MB
- 平均を 16 rank 合計へ換算: **約 81.81 GiB**

これは Python 側の phase ロガーが `after_run_sim` で観測した 74.28 GiB より大きい。
online mode では、NEURON が構築したモデルを同一プロセス内で CoreNEURON 表現へ
転送する。転送・setup 中に両表現や転送バッファが一時的に重なることが、このピーク
の有力な説明である。ただし、保存ログだけから約 81.81 GiB の内部内訳までは確定
できない。

したがって運用上は次のように解釈する。

- 定常的な simulation 終了時 RSS: 約 74.3 GiB
- setup を含む安全側の実測上限: 約 81.8 GiB
- メモリ制限を設定する場合の推奨余裕: 少なくとも 90 GiB、可能なら 100 GiB 以上

今回のジョブは `MemoryHigh=105G`、`MemoryMax=115G` で実行し、上限には到達しなかった。

RSS 合計は共有ライブラリなどの共有ページを process ごとに重複計上し得る。一方、
system available memory は同一マシン上の他プロセスの影響を受ける。そのため絶対的な
物理メモリ量ではなく、同じ16 rank・同じホスト上の相対比較として扱う。

## 4. seed・入力・LC/NA 条件

全比較で以下を固定した。

| 項目 | 値 |
|---|---:|
| connection seed | 4321 |
| stimulation seed | 1234 |
| location seed | 4321 |
| scale | 0.3 |
| `ihGbar` | 0.75 |
| temperature | 34 °C |
| initial Vm | -80 mV |
| integration dt | 0.025 ms |

接続生成は CoreNEURON への切替前に NetPyNE の Random123 で完了し、June NEURON
と CoreNEURON の rank 別接続数も一致した。MOD 内に乱数呼出しはなく、NetStim、
IClamp、pulse は無効である。

3 s 全実行の長距離 VecStim は 12,562 発で、NEURON 9.0.1、NEURON 8.2.7、
CoreNEURON の全てで GID・時刻・順序が完全一致した。matched-prefix の先頭
100 ms でも外部入力 416 発が完全一致した。

明示的な LC population や LC 専用入力はない。以前 LC/NA 相当として扱われた主要
パラメータは HCN 系コンダクタンスをスケールする `ihGbar=0.75` であり、全条件で
同一である。動的 ih/hd 変更も無効である。TPO、TVL、S1、S2、cM1、M2、OC の
長距離入力レートも同一である。

## 5. 段階的な同一性検証

### 5.1 発火なし100 ms条件

同じ100 ms入力生成条件、同じ互換化MODで比較した。

| 指標 | NEURON | CoreNEURON |
|---|---:|---:|
| 外部入力 | 81 | 81 |
| 物理細胞スパイク | 0 | 0 |
| 全PT5B Vm | bitwise equal | bitwise equal |
| Vm RMSE | 0 mV | 0 mV |

この結果から、基本的なモデル転送、固定刻み積分、Vm記録、VecStim直列化は静かな
状態では一致することが確認できた。

### 5.2 3 s入力系列を用いた発火あり先頭100 ms

| 指標 | NEURON | CoreNEURON |
|---|---:|---:|
| 外部入力 | 416 | 416 |
| 物理細胞スパイク | 586 | 16,270 |
| 全スパイク | 1,002 | 16,686 |

外部入力列はビット単位で一致した。互換化済みNEURONの全1,002スパイクと全PT5B
Vmは、June NEURON基準ともビット単位で一致した。

CoreNEURONのみ、最初の物理スパイクから発火GID・時刻が異なる。全PT5B Vmは
47.0 msから差が現れ、差は以下の速度で増幅した。

| 最大PT5B Vm差の閾値 | 最初に超えた時刻 |
|---|---:|
| 1e-6 mV | 47.0 ms |
| 1e-4 mV | 47.2 ms |
| 1e-3 mV | 47.5 ms |
| 0.01 mV | 48.3 ms |
| 0.1 mV | 50.5 ms |
| 1 mV | 57.3 ms |
| 10 mV | 74.0 ms |

これは、長時間積分で小さな丸め誤差がカオス的に増幅しただけではない。最初の
ネットワーク発火・イベント配送段階ですでに異なる挙動である。

## 6. CoreNEURON の機構上の違い

CoreNEURON は NEURON が構築したモデルを別の計算表現へ転送し、同種のチャネル・
シナプスをまとめて並列実行する。NEURON の Array of Structures に対して
CoreNEURON は Structure of Arrays を用い、SIMD/vectorization とメモリ帯域効率を
改善する。

公式資料が示す主要な互換性境界は次の通り。

- 複数の機構インスタンスが SIMD、thread、GPU で並列実行される
- 実行中に書き込む NMODL `GLOBAL` は共有競合を避けるため `RANGE` が必要
- HOC/C/Python 側データを指す `POINTER` は自動転送できず、必要なデータには
  `BBCOREPOINTER` と直列化処理が必要
- CoreNEURON は Random123 のみをサポート
- per-timestep Python/HOC callback はサポートされない
- 古い MOD の明示的 ion 更新や NEURON 内部 API を使う `VERBATIM` は修正が必要

本実装では、書込み可能な速度論変数の `GLOBAL -> RANGE`、VecStim の
`BBCOREPOINTER` 対応、不要な明示的 `cai` 更新の削除、非互換 `VERBATIM` の整理を
行った。これらの同じ互換化済み MOD を通常 NEURON で実行すると June 基準を完全
再現するため、ソース変更自体は観測差の原因ではない。

関連する公式・一次資料:

- [CoreNEURON Compatibility](https://www.neuronsimulator.org/en/8.2.1/coreneuron/compatibility.html)
- [Running a CoreNEURON simulation](https://www.neuronsimulator.org/en/latest/coreneuron/running-a-simulation.html)
- [CoreNeuron Inputs](https://www.neuronsimulator.org/en/latest/coreneuron/inputs.html)
- [Kumbhar et al. 2019, CoreNEURON: An Optimized Compute Engine](https://doi.org/10.3389/fninf.2019.00063)
- [Official NEURON/CoreNEURON comparison tests](https://github.com/neuronsimulator/testcorenrn)

Kumbhar et al. は適合モデルで NEURON との binary result compatibility を報告して
いる。このモデルで観測された 586 対 16,270 という差は、通常期待される
CoreNEURON の数値差ではない。

## 7. 現時点で疑われる箇所

### 強く疑われる

1. 発火開始時に使われる一部の機構状態が CoreNEURON 表現へ正しく転送されていない
2. 外部 presynaptic spike 自体は一致するが、NetCon、weight、delay、target synapse
   の転送またはイベント配送が一部異なる
3. 詳細細胞または抑制性細胞の一部機構に、静かな条件では現れない CoreNEURON 固有
   の非互換性が残る
4. NEURON 8.2.7 wheel に含まれる CoreNEURON バナーが 8.2.2 であることに関連する
   旧バージョンの問題

### 検証により除外または可能性が低い

- seed の違い
- LC/NA 相当パラメータの違い
- 外部 VecStim 時刻列の違い
- 接続生成 seed・rank 数の違い
- `GLOBAL -> RANGE` などの MOD ソース変更そのもの
- GPU、OpenMP thread 数、`cell-permute`
- LFP 計算の影響
- 単純な長時間丸め誤差の増幅

接続数は一致するが、完全な接続辺リストは保存していないため、転送後 NetCon の
全フィールドがバイト単位で同一かは未確認である。

## 8. 速度比較

### 3 s全実行

| 指標 | NEURON | CoreNEURON |
|---|---:|---:|
| 総時間 | 20,718 s | 19,345 s |
| simulation時間 | 約20,480 s | 19,065 s |
| 全スパイク | 135,093 | 225,885 |

見かけ上は総時間が約6.6%短いが、CoreNEURON側は90,792発多くスパイクを生成して
おり、同じ計算結果ではない。この値を純粋なbackend speedupとして扱うべきではない。

### 結果が完全一致した発火なし100 ms条件

| 指標 | NEURON | CoreNEURON | speedup |
|---|---:|---:|---:|
| simulation時間 | 1,654.28 s | 541.50 s | 3.05x |
| 総時間 | 1,877.21 s | 756.67 s | 2.48x |

この条件では出力が完全一致するため、CoreNEURONの計算エンジンとしての速度上の
利点は確認できる。ただし発火・再帰イベントがないため、発火ありネットワークの
正しい性能値ではない。

## 9. 推奨する次の診断

1. 35--55 ms に限定し、最初に発火する全細胞の soma Vm、閾値、主要イオン電流、
   synaptic conductance を両backendで記録する
2. 最初に不一致となるGIDのcell type/populationを特定する
3. cell typeまたは機構を段階的に無効化し、CoreNEURON差が消える最小構成を探す
4. NetConのpreGid、postGid、weight、delay、target mechanismを転送前後でhash比較する
5. NEURON/CoreNEURON 9系のsource buildで同じmatched-prefix試験を実行する
6. 最小再現モデルが得られたら、NEURON公式issueへ報告する

修正が確認されるまでは、発火あり本番計算には通常NEURONを使用し、CoreNEURONは
診断・性能試験用に限定する。

## 10. 関連実装・成果物

- CoreNEURON MOD build: `scripts/build_coreneuron_mechanisms.py`
- full network runner: `scripts/run_full_quiet_pt5b.py`
- backend comparison: `scripts/compare_neuron_backends.py`
- 実装説明: `docs/coreneuron_backend.md`
- CoreNEURON manifest:
  `docs-temp/full-model-runs/full-quiet-3s-coreneuron827-16r-trial0/manifest.json`
- 全体比較JSON:
  `docs-temp/full-model-runs/full-quiet-3s-coreneuron827-16r-trial0/comparison_to_neuron9.json`

関連コミット:

- `5e3bfe9 Add CoreNEURON full-network backend`
- `9624ee6 Skip legacy globals for CoreNEURON range variables`
- `82019b0 Add matched-prefix backend comparison`
