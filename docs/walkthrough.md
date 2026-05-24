# Patterned Optogenetic Stimulation Simulator Walkthrough

## 目的

Thy1-ChR2 マウスの大脳皮質 M1 に対する高速パターン光刺激を、実験前に計算機上で探索できるようにする。最終的には、window 内の 2 um 程度のピクセル ON/OFF パターン列から、L5B を中心とした皮質ネットワーク状態や運動出力に近い指標を予測する。

この計画では、最初から Dura-Bernal et al., 2023 のフルスケールモデルを直接回すのではなく、次のように分解して軽量化する。

1. 単一 L5B 細胞における ChR2 光応答モデルを作る。
2. 光刺激パターンを、各細胞への入力変化または発火確率に変換する。
3. 変換後の入力を、空間情報を落としたネットワークグラフへ注入する。
4. 複数の初期状態に対してロバストに似た応答を作る刺激パターン列を探索する。

## 全体 Walkthrough

### 1. 実験系を計算問題に分解する

実験では、直径 3 mm の glass window 内に内接する正方形領域へパターン光刺激を投影する。2 um 四方のピクセルで考えると、理想的にはおよそ 1500 x 1500 pixel 規模になる。ただし実際の計算でこの全ピクセルをそのまま扱うと重すぎるため、段階的に粗視化する。

最初に分けるべき情報は以下。

- 空間情報が必要な部分: 光がどの樹状突起・細胞領域に当たるか、window 内の細胞位置、光の減衰、照射点と細胞応答の関係。
- 空間情報が不要になる部分: 光入力が細胞ごとの電流、発火確率、状態変化へ変換された後のネットワーク伝播。

この分解により、重い形態シミュレーションは「光刺激から細胞入力への変換関数」を作るために使い、ネットワーク探索はより軽いグラフ・力学系として扱う。

### 2. 単一細胞モデルを作る

最初の実装対象は L5B pyramidal cell 1 個である。Dura-Bernal et al., 2023 の公開モデルから、M1 L5B に相当する形態・チャネル設定を抽出する。もし複数候補がある場合は、代表的な数個を選び、細胞間ばらつきとして扱う。

この単一細胞モデルに ChR2 conductance を追加する。Thy1-ChR2 では soma だけでなく dendrite にも発現している想定なので、初期モデルでは全コンパートメント発現、apical/basal/soma 別発現、dendrite 優位発現などをパラメータとして振る。

単一細胞シミュレーションで得たいものは以下。

- 光の照射位置と発火閾値の関係。
- 光強度、照射時間、PWM duty cycle と膜電位応答の関係。
- 単発発火、burst、failure の条件。
- 光誘発スパイク波形が自然発火波形からどの程度ずれるか。
- Na channel inactivation などにより、高頻度刺激や長時間照射で応答がどう変化するか。

### 3. 細胞ごとの光受容野を推定する

単一細胞モデルに対して、脳表からの照射点を細胞中心からずらしながら応答を測定する。これにより、各細胞について「どの空間位置に光が当たると、どれだけ入力電流・膜電位変化・発火確率が生じるか」という光受容野を作る。

初期近似では、光による ChR2 電流は足し合わせ可能とみなす。つまり、複数ピクセルの同時照射は、各ピクセルからの電流寄与の線形和として扱う。ただし発火・burst・Na channel inactivation は非線形なので、線形化する対象は「光から conductance/current への変換」までに限定する。

出力としては、次のような変換を作る。

```text
pattern pixels -> cell-specific ChR2 conductance/current -> firing probability or state perturbation
```

### 4. ネットワークをグラフ化する

Dura-Bernal et al. のフルモデルから、細胞種、接続、シナプス重み、伝達遅延、時定数、短期可塑性などを抽出する。ここでは細胞の 3D 位置は、光入力を割り当てる段階では使うが、ネットワーク伝播そのものでは落としてよい。

グラフの node は細胞または cell-type population、edge は synaptic influence とする。edge には最低限以下を持たせる。

- pre/post cell type
- excitatory/inhibitory
- synaptic weight
- delay
- rise/decay time constant
- connection probability or multiplicity
- short-term plasticity parameters, if available

このグラフ上で、光刺激による細胞状態変化を外部入力として与える。

### 5. 初期状態にロバストな刺激パターンを探索する

実験上の問題は、同じ刺激でも初期状態により応答が変わることである。したがって最終的な最適化問題は、特定の初期状態で大きな応答を出すことではなく、複数の初期状態に対して似た出力軌道へ収束させる刺激を探すことになる。

候補となる評価関数は以下。

- 複数初期状態間の出力ばらつきが小さい。
- 前肢・後肢関連 population の活動が目的の位相差で交互に変化する。
- 過剰同期や異常 burst を避ける。
- 総光量、最大光量、連続照射時間を制約内に保つ。
- 実験装置の 9000 Hz 切り替え制約を満たす。

## 単一細胞シミュレーション Plan

### Phase 0: 環境と前提の固定

Python を使う場合は、プロジェクト内で `uv venv --python=3.13` により仮想環境を作成する。NEURON、解析用ライブラリ、可視化ライブラリはこの仮想環境内に入れる。

この phase で決めること。

- NEURON Python API を使う。
- 形態ファイル、MOD ファイル、パラメータファイルの配置規則を決める。
- Dura-Bernal et al. モデルから抽出する L5B 細胞候補を決める。
- 実験で使う光波長、最大光強度、パルス幅、PWM 周波数の候補範囲を決める。
- 単位系を固定する: um, ms, mV, nA, mW/mm2 など。

### Phase 1: ベースライン L5B 細胞を再現する

ChR2 を入れる前に、L5B 細胞単体として妥当な挙動を確認する。

確認項目。

- resting membrane potential
- input resistance
- membrane time constant
- rheobase
- current injection に対する f-I curve
- single spike waveform: amplitude, half-width, AHP
- burst tendency

この段階で元モデルの再現性が取れない場合、ChR2 応答の解釈ができないため先に修正する。

### Phase 2: ChR2 mechanism を追加する

ChR2 を NEURON mechanism として追加する。最初は実験に完全一致した詳細モデルよりも、パラメータ探索しやすい conductance-based model を優先する。

最小モデルの形。

```text
I_ChR2 = g_ChR2(t, irradiance) * (V - E_ChR2)
```

必要なパラメータ。

- reversal potential
- maximal conductance density
- activation time constant
- deactivation time constant
- desensitization or inactivation term
- recovery time constant
- irradiance-to-open-probability conversion

発現分布の候補。

- soma + dendrite uniform
- dendrite-only
- apical-dominant
- basal-dominant
- distance-from-soma dependent
- compartment-type dependent

### Phase 3: 光伝播モデルを最小実装する

最初の目的は正確な光学シミュレーションではなく、位置依存の入力を細胞へ与えられるようにすることである。

初期モデルでは、脳表上の照射点から各 compartment への寄与を以下で近似する。

```text
irradiance(compartment, t) =
    surface_pattern(x, y, t)
    * lateral_spread(distance_xy)
    * depth_decay(depth_z)
```

候補。

- lateral_spread: Gaussian
- depth_decay: exponential
- pixel integration: compartment midpoint で近似し、必要なら segment 長で重み付け

この phase では、2 um pixel をそのまま全探索せず、まず 20-50 um 間隔の粗い照射点 grid で応答を調べる。

### Phase 4: 空間応答マップを作る

細胞中心を基準に、照射点をずらして応答を測る。

掃引軸。

- x-y offset from response center
- spot size
- irradiance
- pulse duration
- compartment expression pattern
- cell morphology candidate

出力。

- subthreshold peak depolarization
- spike probability
- first spike latency
- number of spikes
- burst probability
- charge or integrated ChR2 current

この結果から、単一細胞の光受容野を作る。

### Phase 5: 時間応答と PWM 応答を調べる

照射位置を固定し、時間パターンを掃引する。

掃引軸。

- single pulse duration
- pulse train frequency
- PWM carrier frequency
- PWM duty cycle
- sequence duration
- inter-stimulus interval
- total light dose

解析項目。

- 発火閾値となる最小 pulse duration
- 発火 latency のばらつき
- burst への移行条件
- spike waveform distortion
- Na channel inactivation の蓄積
- ChR2 desensitization の蓄積
- 同じ総光量で continuous と PWM がどう違うか

### Phase 6: 実験データに合わせてパラメータを推定する

in vivo patch clamp または先行実験から得られるデータに合わせ、ChR2 発現量・光伝播・kinetics を調整する。

合わせる観測量。

- 光強度に対する depolarization
- 光強度に対する spike probability
- pulse duration threshold
- latency distribution
- spike count distribution
- voltage trace shape
- repeated stimulation での adaptation

推定対象。

- ChR2 maximal conductance density
- expression distribution
- optical spread
- depth attenuation
- ChR2 activation/deactivation/desensitization

最初は global optimization ではなく、grid search と Bayesian optimization の中間程度の現実的な探索でよい。

### Phase 7: ネットワーク入力用の surrogate を作る

NEURON の単一細胞モデルを毎回ネットワーク最適化に使うと重いので、入力変換用の surrogate を作る。

候補。

- lookup table
- generalized linear model
- small neural network
- Gaussian process
- low-rank receptive field + temporal filter

最初の実装では lookup table を推奨する。理由は、物理量との対応が明確で、後から in vivo patch clamp データで修正しやすいからである。

## 最初に作るべき最小成果物

レビュー後に実装へ進む場合、最初の milestone は以下にする。

1. `uv venv --python=3.13` で Python 環境を作る。
2. NEURON Python が動くことを確認する。
3. 単純な soma-only または ball-and-stick cell で ChR2 mechanism を動かす。
4. 1 pulse 光刺激で depolarization と spike が起きることを確認する。
5. pulse duration と irradiance の 2D sweep を行い、発火閾値マップを出す。
6. Dura-Bernal 由来 L5B morphology へ置き換える。
7. 空間 offset sweep を追加し、光受容野の初版を作る。

この順序にすると、Dura-Bernal モデル抽出で詰まっても ChR2 と光刺激部分の開発は前に進められる。

## 現時点の設計上の注意点

- 「光入力から ChR2 電流」は線形近似できる可能性が高いが、「電流から発火」は明確に非線形である。
- 9000 Hz のパターン更新をそのまま NEURON に入れると計算が重くなるため、PWM は有効 conductance 近似と明示的 pulse train の両方で比較する。
- 麻酔下 M1 の運動誘発は初期状態依存性が強いので、単一細胞モデルだけで運動を説明しようとせず、単一細胞はあくまで光入力変換器として位置づける。
- Thy1-ChR2 の発現量・分布は不確実性が大きいため、固定値ではなく推定パラメータとして扱う。
- Dura-Bernal フルモデルの詳細な biophysics を残しすぎると重くなるため、単一細胞で必要な機能と、ネットワークで必要な機能を分離して実装する。

