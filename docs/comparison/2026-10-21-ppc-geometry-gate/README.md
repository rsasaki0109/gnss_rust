# PPC FLOATのGDOP・共分散・原観測診断

2026-10-21。前回の最大約938 mの外れ値について、幾何の弱さを調べた。
GDOP 50の追加ゲートを評価したが、最大外れは解消しなかった。
ゲートは既定で無効のまま。実用性能の受け入れ結果ではない。

## 実装と独立検証

geometry::dilution_of_precisionは単一共通時計のH=[-LOS,1]に対する
無重みのGDOP/PDOP/TDOPと幾何共分散を返す。
既存のHouseholder QRを共有し、正規行列の逆算や対角の下限補正は行わない。
これはRTKの事後共分散とは別の診断値で、位置精度の保証ではない。

固定した元C++ ObservationData::calculateDOPを呼ぶ16ケースと、
独立のEigen逆行列による幾何共分散を比較した。
4/5/7/9衛星、北/南半球、通常/弱い幾何を含む。
通常ケースは相対1e-10、弱い幾何は元の正規行列逆算の精度を考慮して相対1e-6。
行順の変更、重複衛星・ランク不足・非有限入力も検証した。

RtkFloatConfig::max_geometry_gdopを追加し、Noneなら従来の採否を維持する。
有効な場合は最終推定位置のGDOPを確定前に検証し、超過をGeometryTooWeakとして棄却する。
数値状態は確定せず、次の採用時は既存の失敗後ルールに従ってアークを初期化する。
上限との厳密な一致は通過し、1 ULP下の上限では棄却する境界テストを含む。

## 凍結した比較

全6走行に対してゲートなし／GDOP 50の2設定を各1実行、計12実行した。
50は固定したSPPProcessor::SPPConfigの既定値から採った実験値で、
元RTK政策の同一性や、PPCの真値から決めた受け入れ閾値ではない。
ソース・入力・設定・バイナリを [predeclared-inputs.json](predeclared-inputs.json) に
保存してから実行を開始した。

ゲートなしの全6実行で、前回e6502f4のCSVの既存18列がすべて文字列として一致した。
新しくGDOP、PDOP、位置共分散、使用衛星の診断列だけを加えている。
今回もGPS 1Cの1 Hz同時刻対だけを使い、片側だけのエポックやINVALIDを保持した。
FIX・PPP・元C++全RTKソルバーは実行していない。速度・メモリの受け入れ評価もない。

| 走行 | ゲートなしFLOAT | GDOP 50 FLOAT | なしの水平RMS m | 追加後の水平RMS m | 追加後の最大3D誤差 m |
| --- | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 840 | 840 | 2.102 | 2.102 | 99.155 |
| nagoya-run2 | 1241 | 1221 | 25.083 | 24.885 | 937.912 |
| nagoya-run3 | 697 | 613 | 13.092 | 12.234 | 109.199 |
| tokyo-run1 | 1604 | 1451 | 11.109 | 8.854 | 207.093 |
| tokyo-run2 | 1471 | 1467 | 4.735 | 4.732 | 48.475 |
| tokyo-run3 | 2405 | 2392 | 9.868 | 9.849 | 189.910 |

FLOAT出力は全11,733対中8,258→7,984（70.4%→68.0%）。
位置誤差の表はreference.csvをアンテナ位置とみなす仮定で、各設定の全FLOAT出力を使用。
両設定で出力する時刻が異なるため、RMSの差をそのまま精度改善と判断できない。
IMU位置・Tokyoの両レバーアームを含む全仮定も [report.json](report.json) に保存した。
ゲートを出力後の真値フィルターとして使う処理はない。

## 最大外れの調査

[worst-positions.json](worst-positions.json)に各設定の上位10外れを、
時刻・使用衛星・GDOP・共分散・残差とともに保存した。
これらの抽出は実行後の調査専用で、ソルバーや採否に戻していない。

Nagoya run2のGPS週2323/週内秒556760は、ゲートなしで3D誤差937.912 m、
GDOP約10.009、事後位置共分散のトレース平方根約4.221 mだった。
Tokyo run1の208.317 mの外れでもGDOPは約6.924だった。
最大級の外れは高いGDOPに限定されず、推定された共分散だけでも説明できない。

両点の前後を原観測で調べる読込専用例も追加した。
コード−位相の単差、前回からの変化、Dopplerとの位相差分誤差、LLI、C/N0を保存する。
調査窓は結果から選んだもので、事前に決めた性能比較とは分けている。
結果と出典は [raw-phase-probes.json](raw-phase-probes.json) とCSVのgzip。
衛星ごとのambiguity定数があるため、CMCの絶対値だけで異常とは判定していない。

最悪点では、新しく共通の位相が利用可能になったG32の移動局C/N0が25.938 dB-Hz。
Tokyoの外れでも新しいG07は29.313 dB-Hzだった。
いずれもその時刻には比較できる前回のCMC/Doppler履歴がなく、初回利用の品質確認が
必要なことが分かった。低C/N0が原因と確定した結果ではなく、
今後は再捕捉時の観測採用・ambiguity初期化と継続性を調べる。
SNRやCMCの新しい閾値は今回の測位実行には適用していない。

## 保存と再検証

生のCSV/標準エラーをgzipで保存し、実行コマンド・終了コード・非圧縮SHA-256を記録した。
元C++の16ケースは再生成して完全一致を確認。
追加4テストを含むRust276（docを含む）、examples、fmt、全ターゲットclippy、
releaseビルドが通過した。[qualification.json](qualification.json) に記録した。

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_float
python3 -B scripts/run_ppc_geometry_gate.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float \
  --output-dir /tmp/new-ppc-geometry-gate
~~~

保存した証拠の検証は次で行う。--dataset-rootを付けると参考誤差も再計算する。
ソルバーやネットワークは実行しない。

~~~bash
python3 -B docs/comparison/2026-10-21-ppc-geometry-gate/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

PPCでのRTK・PPP同等性能は未達で、確認済みチェックリストは18/24（75%）。
