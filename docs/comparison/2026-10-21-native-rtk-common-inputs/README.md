# PPCでの元C++全RTKソルバーとRustの初回比較

2026-10-21。固定したlibgnss++のRTKProcessorを新しくビルドし、PPC全6走行を
Rustと同じGPS 1Cの同時刻対で最後まで実行した。
**これは制限した設定での診断比較で、PPCでの同等性能の受け入れではない。**
元C++にも大きな外れとSPPへの退避があり、比較元の実用設定はまだ確立できていない。

## 入力と設定

元ソースは72f3b7c2c4c088dc575286a6f0edf4e407bc499a。
チェック済みアーカイブに対して2,091ファイルを照合し、変更せずReleaseビルドした。
[native-build.json](native-build.json)にコマンド、コンパイラ、実際のflagsを保存。
CMakeのPython起動スクリプトはimportに失敗したため、同梱の実バイナリで再実行した。
失敗ログと成功したconfigure/buildログを両方残している。

元C++のRINEXReaderとRustの読込器から、GPS 1Cのコード・位相・Doppler・C/N0・LLIを
別々に出力した。全6走行の184,453行について、数値のf64ビット表現、コード/位相/
Dopplerの有無、未対応エポックから引き継ぐreset通知が一致した。
元Observation APIに存在しないSNR/LLIの欠測ビットだけは0へ正規化している。
これは**観測値の照合**であり、航法データのデコード値やソルバー内部の採用衛星・
補正モデルが一致したことを示すものではない。同じbase.navを入力した。

GPS週を含む時刻差1e-7秒以内を対応付け、片側だけのエポックと両末尾も全て出力。
全11,733対、移動局だけ46,901、基準局だけ13。元C++の有効出力に古い時刻の解はなかった。
PPCの原エポックflagは全て0。一般のRINEXイベントに対応した比較器とはしていない。

両者ともKINEMATIC、仰角10度、基線上限10 km、測量済み基準局位置を使う。
最初は移動局ヘッダー位置、以後は自分自身の直前の有効RTK位置を与える。
元C++ではこのために**既定falseのprefer_rover_position_seedをtrue**にした。
元C++の通常のSPPによる初期化を使う設定とは異なり、次に切り分ける対象である。
元C++のARは既定CONTINUOUS/EXTENDEDのまま、QC・スリップ検出・退避も有効。
Rustは現状GPS単周波FLOATのみで、移動用FIX政策はない。
同じ内部処理を実装したという比較ではない。

未対応エポックのLLI loss/half-cycle通知は次の対応対で全アークを再初期化する。
元C++へは観測のコピーのLLI bit 1に反映し、照合台帳には変更前の原LLIと通知列を残す。
真値をソルバー・初期値・採否判定へ渡す処理はない。
その他の差を含む凍結設定とソース/入力/バイナリ/4静的ライブラリのSHA-256は
[predeclared-inputs.json](predeclared-inputs.json)に実行前に保存した。

## 結果

| 走行 | 対応対 | 元C++ FLOAT | FIX | SPP退避 | INVALID | Rust FLOAT |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 1521 | 350 | 3 | 519 | 649 | 840 |
| nagoya-run2 | 1891 | 219 | 3 | 604 | 1065 | 1241 |
| nagoya-run3 | 1041 | 70 | 16 | 347 | 608 | 697 |
| tokyo-run1 | 2388 | 273 | 0 | 1070 | 1045 | 1604 |
| tokyo-run2 | 1831 | 161 | 20 | 1143 | 507 | 1471 |
| tokyo-run3 | 3061 | 242 | 2 | 2082 | 735 | 2405 |
| 合計 | 11733 | 1315 | 44 | 5765 | 4609 | 8258 |

元C++のRTK出力は1,359/11,733（11.6%）、Rustは8,258/11,733（70.4%）。
SPPはRTK成功へ算入していない。元の5 Hz移動局58,634エポックのうち対応対だけを
処理しているため、この割合を全入力に対する継続率とは扱わない。
Rustの全6 CSVは前回b711684のゲートなし出力と全列が一致した。

| 走行 | 元C++ RTK水平RMS m | Rust水平RMS m | 元C++最大3D m | Rust最大3D m |
| --- | ---: | ---: | ---: | ---: |
| nagoya-run1 | 74.363 | 2.102 | 147.966 | 99.155 |
| nagoya-run2 | 79.781 | 25.083 | 965.496 | 937.912 |
| nagoya-run3 | 181.405 | 13.092 | 584.067 | 109.137 |
| tokyo-run1 | 1180.363 | 11.109 | 3132.849 | 208.317 |
| tokyo-run2 | 47.694 | 4.735 | 146.697 | 48.475 |
| tokyo-run3 | 46.371 | 9.868 | 149.179 | 189.910 |

表はreference.csvをアンテナ位置とみなす参考評価。各ソルバーの全RTK出力を使い、
成功時刻が異なるのでRMSだけで優劣を判断しない。
IMU仮定とTokyoの異なる二つのレバーアームも全て[report.json](report.json)に保存。
[diagnostics.json](diagnostics.json)には共通成功時刻/片側成功/両側不成立の数、
共通時刻だけの誤差、元C++のFLOAT/FIX/SPP別誤差、スリップ検出数を保存した。
例えばNagoya run1の3 FIXはこの仮定で水平RMS約124 m、Tokyo run3の2 FIXは約50 m。
FIXという状態だけでは正しさを保証できない。誤FIXの受け入れ閾値はまだ定めていない。

元C++の設定・単周波制限・初期位置更新の影響を先に調べる必要がある。
この結果で「比較元より高性能」「同等性を達成」とは主張しない。
全衛星系/二周波、移動用FIX、PPP製品の信号基準への接続、反復した速度/メモリ評価も残る。

## 検証と再実行

移動合成64エポックのsmokeでは元C++が45 FIX・19 FLOAT、最大3D誤差約0.84 mm。
1,152入力行もRustと一致。[synthetic-smoke.json](synthetic-smoke.json)とCSVを保存した。
これはラッパーの基本動作の確認であり、都市環境での精度の保証ではない。
Rust276テスト、examples、fmt、全ターゲットclippy、releaseが通過。
[qualification.json](qualification.json)にコマンド・終了コード・圧縮ログのハッシュを記録。

最初の実行は、対応しない小数時刻の異なる十進表記をDecimalで比較するガードが失敗した。
両表記は同一f64だったため、f64ビット表現の比較に修正した。
ソルバー/閾値の変更はなく、[失敗した実行](../2026-10-21-native-rtk-common-inputs-format-attempt/README.md)
も別ディレクトリに保存して全6走行を新たに実行した。

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_float --example ppc_input_ledger
# native-build.jsonのコマンドで元C++とadapterを構築した後:
python3 -B scripts/run_native_ppc_comparison.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --native-binary /workspace/.cloud-setup/solver-reference-20261021/native_ppc_rtk \
  --native-build /workspace/.cloud-setup/solver-reference-20261021 \
  --rust-binary target/release/examples/ppc_rtk_float \
  --rust-ledger-binary target/release/examples/ppc_input_ledger \
  --output-dir /tmp/new-native-ppc-comparison

python3 -B docs/comparison/2026-10-21-native-rtk-common-inputs/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

CSV・stderr・終了コード・ハッシュは保存済み。nativeのCSVは明示したファイルから取得し、
プロセスstdoutは別ファイルに保持した。観測台帳は初回に行数・正規化SHA・両原CSVのSHAを
保存した後、同じバイナリで再実行して両台帳もgzipで保持した。
再実行した元C++の全CSVと両台帳のSHAは初回に完全一致した。
verifyは保持した両台帳の184,453行を独立に比較し直す。
再実行のコマンドと終了コードは各input-retention-replay.jsonに保存した。
観測台帳の出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)。
Copyright (c) 2025 Taro Suzuki、[MIT License](LICENSE-PPC-Dataset)を同梱した。
入力のZIP・各ファイルのハッシュと取得元は
[入力検証記録](../2026-10-21-ppc-input-readiness/README.md)を参照。
将来ソースが変わった場合、この証拠を追加したコミットのソースでverifyを実行する。

PPCでの実用RTK・PPP同等性能は未達。確認済みチェックリストは18/24（75%）。
