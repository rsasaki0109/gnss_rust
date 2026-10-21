# 初期残差3 mのエポック保留ゲート

2026-10-21。初回/再捕捉衛星にはCMC/Doppler履歴がないため、初期残差も評価した。
元RTKConfig::outlier_threshold=3 mを使う選択式ゲートを追加し、全6走行・6設定を実行。
**出力率が大きく落ち、SPPと併用しても最大約1,573 mの外れは残った。既定では無効。**
同等性能や精度改善の受け入れ結果ではない。

## 元実装との関係

元C++のsuppressOutlierRowsは閾値を超えた行のH/vを0にし、相関したRは変更しない。
これは、その行を削除してRの対応する主部分行列を使う処理とは異なる。
1状態・2行・Rの相関0.5の閉形式反例で、残した測定1に対し、前者の事後平均/分散は
4/7と3/7、後者は1/2と1/2になることをテストした。

今回のRust実装は行の分類だけを照合し、**超過時はエポック全体の更新を保留する**。
H/v/Rや閾値を変更して測定を通す処理はない。元C++全体の外れ値政策の移植ではない。
元suppressOutlierRowsを直接呼ぶ10ケースで、符号・3 mとの一致・1 ULP内外とマスクを照合。
有効入力の行分類は完全一致し、無効入力はRustで明示的に拒否する。

RtkFloatConfig::max_prefit_residual_m=Some(3.0)で、SPP/従来位置による初期化と
ambiguityリセットの後、最初の線形化時に1回だけDDコード/位相残差を検査する。
どれかが厳密に3 mを超えるとPrefitResidualRejectedを返し、数値状態・履歴を確定しない。
次の採用時は従来の失敗後の全アーク再初期化に従う。
last_prefit_diagnosticは成功/拒否の両方で初期値の出典・座標、最大残差、拒否したDD対、
当該試行のスリップ証拠を保持する。新しい試行ではまず消去し、以前の診断を混同しない。
真値をゲートや設定選択に渡す処理はない。

追加11テストで元10分類、相関Rの反例、初回衛星の200 mコード外れ、位相ジャンプ、
状態不変・復帰、閾値の一致/1 ULP、SPP初期化との順序、無効設定を検証した。
非制限的なゲートを使うと全32合成エポックの数値状態は従来と完全一致する。

## 全6走行・36実行

仰角・ノイズ・基線上限・NISなど他の閾値は変更せず、次の6設定を凍結した。
GPS 1Cの11,733対応対と全片側入力を出力し、元5 Hz全時刻の測位とは扱わない。

| 走行 | 従来FLOAT | SPPのみ | ゲートのみ | SPP＋ゲート | 検出＋ゲート | SPP＋検出＋ゲート |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 840 | 840 | 69 | 558 | 69 | 557 |
| nagoya-run2 | 1241 | 1241 | 55 | 795 | 55 | 797 |
| nagoya-run3 | 697 | 697 | 38 | 281 | 38 | 282 |
| tokyo-run1 | 1604 | 1599 | 20 | 806 | 39 | 806 |
| tokyo-run2 | 1471 | 1471 | 45 | 804 | 45 | 809 |
| tokyo-run3 | 2405 | 2405 | 41 | 1461 | 41 | 1464 |
| 合計 | 8258 | 8253 | 268 | 4705 | 287 | 4715 |

従来とSPPのみの全CSVは先行記録と全列が一致した。
ゲートだけでは古いRTK位置を初期値に使い続けるため、移動に追従できず出力を大きく失う。
この268点だけの小さい誤差を改善と主張しない。

| 走行 | SPP＋ゲート 水平RMS m | 鉛直RMS m | 最大3D m |
| --- | ---: | ---: | ---: |
| nagoya-run1 | 0.928 | 2.802 | 28.870 |
| nagoya-run2 | 21.174 | 109.696 | 1572.517 |
| nagoya-run3 | 7.865 | 24.959 | 149.374 |
| tokyo-run1 | 8.403 | 21.297 | 236.450 |
| tokyo-run2 | 2.491 | 3.686 | 27.058 |
| tokyo-run3 | 4.763 | 9.292 | 192.539 |

表はreference.csvをアンテナとみなす事後参考評価。各設定の全出力を使い、成功時刻が違う。
IMU仮定とTokyo両レバーアーム、全無効理由・ゲートに到達した失敗試行の診断数も
[report.json](report.json)へ保存した。全体の精度・継続率改善とは判断しない。

Nagoya run2の最大外れは、初期SPP自体が参考座標から約1,572.954 m離れ、
初期残差3 mを通過した後にRTKが約1,572.517 mの誤差で出力していた。
[生出力を含む調査記録](worst-spp-gated-nagoya2.json)に上位10点を保存。
これは事後調査であり、SPPやRTKの採否には戻していない。
小さい残差だけでは、自己整合している誤った初期位置を識別できない。
次はSPP側の衛星数・幾何・外れ値処理と、衛星系を増やす経路を調べる。

## 再実行

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_float
python3 -B scripts/run_ppc_prefit.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float \
  --output-dir /tmp/new-ppc-prefit-gate

python3 -B docs/comparison/2026-10-21-prefit-epoch-gate/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

生CSV・stderr・終了コード・コマンド・SHAを保存。Rust302テスト、examples、fmt、
全ターゲットclippy、releaseが通過した。元C++全RTKの再実行や速度/メモリ受け入れはない。
PPCでのRTK・PPP同等性能は未達。確認済みチェックリストは18/24（75%）。
