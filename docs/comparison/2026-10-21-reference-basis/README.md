# RTK基準衛星変更の診断（2026-10-21）

基準衛星を変えると、二重差の単純RMSは変わる。一方、H・残差・全共分散を
一緒に変換すれば、Gaussian更新の推定値・共分散・innovation NISは変わらない。
pinned C++の22ケースでこの違いを確認した。

## 追加したAPI

rtk_reference.DoubleDifferenceBasisは衛星の実IDと入力順を保持し、
reference-minus-targetの各blockについて新DD = M旧DDを作る。
rebaseは H'=MH、v'=Mv、R'=MRMᵀ を同時に計算する。
phase/code間を含む全cross covarianceを保持し、同じ衛星系の完全なDD blockを扱う。
未知reference、重複/混在ID、不正寸法、非有限、非正定値を拒否する。
これは測定basisの変換であり、DD状態を持つfilterのstate/prior/Qの変換は別に必要。

nominal_weighted_sum_of_squaresは vᵀR⁻¹v。
postfit残差へ使った値はinnovation NISや校正済みchi-squareとは扱わない。

## 検証した二つの違い

3衛星の旧DD phase残差 [0.010, 0.015] mは、別referenceで [-0.010, 0.005] mになる。
RMSは0.012748→0.007906m。既存phase RMS gateの0.01mを跨ぐが、
全共分散での重み付き残差量は同じ。
これは一つのgateの性質を示す。FIXの全条件を満たすという証明ではない。

独立DD process noiseのIは、reference変更後に [[1,1],[1,2]] となりIではない。
一方、均等な独立SD noiseに由来するDD covariance [[2,1],[1,2]] は、
対応するreference変更でも同じになる。
現コピーの静止RTKはSD状態でQを加算している。
前回のaligned RTKで使ったDD Qとの実装照合には、そのコードの復旧が必要。

## 比較元と検証範囲

libgnss++ commit 72f3b7c2c4c088dc575286a6f0edf4e407bc499a の
unchanged measurement assembly／Kalman componentを呼び、4/5/6/7衛星の全referenceを計算。
合成SD入力と、宣言したshared phase/code sourceを使った。
H・v・全R・全posterior mean/covariance・NIS・名目残差量を元の1e-10以下の条件で照合。
phase/code cross covarianceを落とすと残差量が変わる反例、任意の衛星順、
往復、入力不変・境界も検証した。許容差やsolver gateは変更していない。

Rust **250 tests**（追加7）、examples **0 tests**、fmt/clippy/release成功。
examplesの0件は機能検証の合格数には含めない。
[元fixtureと来歴](../../../tests/fixtures/upstream_rtk_reference_basis.csv.json)、
[実コマンドとreturncode](qualification.json)、[集計とsource hash](summary.json)。
実行ログは元のbytesをgzip圧縮して保存した。
新PPC実行0、新native whole実行0。実データの精度改善は主張しない。

## 作業コピーの復旧状況

この実行環境はGitHub mainと同じ ccc9b01 で、前回のCMC実装・比較archiveは見つからない。
local branches/reflog、全remote refs、Gitのunreachable objects、workspace保存領域を確認した。
unreachable blobは383/10byteの2個で、最新コードを含まない。
前回の完成度83%は、このコピーでは再検証できない。
現在のGOALSチェックリストは18/24（75%）。実用RTK・PPP性能目標は未達。

前回の成果を含むバックアップ／別の環境snapshotの保存先を確認する。
今回の診断部品は復旧先へ移せるよう独立moduleにし、既存RTK/PPPの既定動作を維持した。
