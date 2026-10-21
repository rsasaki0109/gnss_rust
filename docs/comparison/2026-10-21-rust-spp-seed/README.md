# RustのSPP初期化接続とPPCでの結果

2026-10-21。元C++では[既定SPP初期化](../2026-10-21-native-default-seed/README.md)が
大きく影響したため、Rustにも明示的な選択肢として接続し、全6走行で評価した。
**RustのPPC外れ値は改善しなかった。既定では無効のままにする。**

## 実装と検証

RtkFloatConfig::kinematic_spp_seedへSome(SppConfig)を指定すると、GPS 1Cコードだけで
そのエポックのSPPを解き、キネマティックRTKの位置平均と衛星選択の初期値に使う。
位置分散900 m²と位置/ambiguityの交差共分散を従来どおり初期化し、継続しているSD
ambiguity平均と共分散は保持する。SPPの共分散を重ねて測定更新へ入れる処理はない。
SPP失敗時は直前のRTK位置/初期位置へ戻り、解に初期値の出典・座標・失敗理由を残す。
静止モードとの組合せと無効なSPP設定は拒否する。
SPP出力をRTKの成功やFIXに置き換える処理はない。

既定Noneの経路は数値演算を維持し、今回も前回2972bf5の全6 CSV・全列が完全一致した。
examples/ppc_rtk_floatのkinematic-sppだけが初期値の診断5列を追加する。
既存モードのCSV形式は維持している。

追加4テストで、元C++から独立生成した移動64エポックの真値への追従、
連続ambiguityの保持、SPP失敗時の数値状態の完全一致、GPS 1C以外の除外、
RTK更新失敗時の状態不変とアーク再初期化、設定拒否を検証した。
合成SPPは生成元と同じSaastamoinen有効・電離層なしで検証し、既存の3 mm許容差を維持した。
都市データの精度や、元C++のSPP/QCアルゴリズム全体の同等性を示すテストではない。

## 凍結した全6走行・12実行

従来のキネマティック初期値と、SppConfig::default()によるSPP初期値を各6実行した。
閾値・仰角・基線上限・ノイズ・時刻対応付け・原観測は変えていない。
今回もGPS 1Cの11,733対応対のみを処理し、片側だけの入力とINVALIDを出力した。
真値は事後評価にだけ使用。ソース・バイナリ・入力のSHAは実行前に凍結した。
元C++のSPP初期化結果は、先行する実行の出力を参照している。

| 走行 | 従来FLOAT | SPP初期値FLOAT | 従来水平RMS m | SPP水平RMS m | SPP最大3D m |
| --- | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 840 | 840 | 2.102 | 2.309 | 141.252 |
| nagoya-run2 | 1241 | 1241 | 25.083 | 27.916 | 1572.517 |
| nagoya-run3 | 697 | 697 | 13.092 | 13.293 | 149.374 |
| tokyo-run1 | 1604 | 1599 | 11.109 | 28.752 | 899.603 |
| tokyo-run2 | 1471 | 1471 | 4.735 | 4.640 | 48.599 |
| tokyo-run3 | 2405 | 2405 | 9.868 | 9.650 | 192.539 |

FLOATは8,258→8,253。SPP設定の採用解のうち8,252エポックでSPPが成功し、
Tokyo run1の1エポックだけNonConvergentによる前回RTK位置へのfallbackだった。
失敗したRTKエポックでは初期値の診断を出力していないため、この数は全SPP試行数ではない。
従来の約938 mの最大外れは、SPP初期値では約1,573 mに悪化した。
従来方式を維持し、SPPの有効化を精度対策として推奨しない。

表はreference.csvをアンテナ位置とみなす参考誤差。各設定の全FLOAT出力を使い、
成功時刻が違うことにも注意する。IMU仮定とTokyoの二つのレバーアームを含む
全仮定と無効理由を[report.json](report.json)に保存した。
元C++もGPS単周波制限下で最大約768 mの外れがあり、実用上の比較基準の確立は継続する。
今回の結果は、初期位置設定だけではRustの大きな外れを解消できなかったことを示す。
新しいCMC/SNR閾値や真値を用いた採否は導入していない。

## 再実行

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_float
python3 -B scripts/run_ppc_spp_seed.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float \
  --output-dir /tmp/new-rust-spp-seed

python3 -B docs/comparison/2026-10-21-rust-spp-seed/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

生CSV・stderrと終了コード・SHAを保存。Rust280テスト、examples、fmt、
全ターゲットclippy、releaseが通過し、qualification.jsonと圧縮ログに記録した。
過去の検証スクリプトが旧Rustソースのハッシュを要求する場合は、当該証拠を追加した
コミットで検証する。今回のverifyは今回のソースと全12実行を確認する。

次は元C++の観測品質処理・未フラグスリップと新規ambiguity初期化を調べる。
移動FIX・全衛星系/二周波・PPP補正基準・速度とメモリの受け入れも未完了。
PPCでのRTK・PPP同等性能は未達。確認済みチェックリストは18/24（75%）。
