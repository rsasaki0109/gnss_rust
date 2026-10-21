# コード−位相差・Doppler検出の移植とPPC評価

2026-10-21。元C++の単周波スリップ検出式を選択式でRTKへ接続した。
全6走行・4設定を実行したが、**最大約938 mの外れは残った。既定設定は変えない。**
PPCでのRTK・PPP同等性能の受け入れではない。

## 実装と独立検証

rtk_slipはCMC=(P_rover−P_base)−λ(L_rover−L_base)の時間変化と、
単差位相距離の変化−{−λ(D_rover−D_base)}×dtを評価する。
元C++のrtk_slip_detection.hppを変更せずに呼ぶ28ケースと照合した。
符号、実測ゼロ、共通コード偏りの相殺、5 m/0.2 mの境界と1 ULP内外を含む。
数値は相対1e-12、検出boolは完全一致。固定元ソースは
72f3b7c2c4c088dc575286a6f0edf4e407bc499a。
元関数が無効入力に返す0/falseはRustでは明示的なErrorにし、欠測をゼロへ変えない。
非有限・非正の波長/閾値/時間・有限入力からの演算overflowも検証した。

RtkFloatConfig::slip_detection=Some(RtkSlipConfig)で、採用予定のGPS 1Cを検査し、
検出した衛星のSD ambiguityだけを再初期化する。その時刻の観測値は捨てない。
RtkSlipDiagnosticに衛星ごとの符号付き残差・判定・Dopplerの有無とエラーを残す。
初回・既知LLI・欠測/再捕捉・前回の更新失敗・長い位相間隔を跨いだ比較は行わない。
失敗エポックは履歴へ入れず、次の採用時は既存の全アーク再初期化規則に従う。

既定Noneでは従来の数値経路を維持する。検出器の5 m・0.20 mとDoppler最大間隔5 sは
元RTKConfigと呼出部から採った。履歴は採用エポックだけを反映するため、
元RTKProcessor全体の履歴更新政策を再現したとの主張はしない。
Dopplerは前回と現在に有効な両局値がある場合に使う。実測0は有効、欠測は未使用、
非有限値はDopplerの診断エラーとして保持し、コード/位相の処理まで捨てない。
5秒超ではDoppler比較を省略して現在の有効値から履歴を再開する。

追加11テストで成分28ケース、LLIなし50/2 cycleジャンプ、対象衛星だけの初期化、
欠測・非有限・ゼロの区別、再捕捉、失敗時の状態不変、5秒境界と復帰、
既知イベント・設定不正を検証した。全スリップを検出できるとの保証ではない。

## 全6走行・24実行

検出なし／CMCのみ／Dopplerのみ／両方を各6実行した。
SPP初期化、追加GDOPゲート、移動FIXは使わず、他の閾値を変更していない。
ソース・バイナリ・入力・設定を実行前に凍結し、真値は事後評価だけに使った。
GPS 1Cの11,733同時刻対と全片側入力を残し、元5 Hz全時刻を測位したとは扱わない。

| 走行 | なしFLOAT | CMC FLOAT | Doppler/両方 FLOAT | なし水平RMS m | 両方水平RMS m | 両方鉛直RMS m | 両方最大3D m |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 840 | 840 | 840 | 2.102 | 2.093 | 7.314 | 99.155 |
| nagoya-run2 | 1241 | 1241 | 1242 | 25.083 | 25.372 | 120.773 | 937.912 |
| nagoya-run3 | 697 | 697 | 700 | 13.092 | 12.599 | 20.001 | 109.117 |
| tokyo-run1 | 1604 | 1604 | 1605 | 11.109 | 11.558 | 22.170 | 208.314 |
| tokyo-run2 | 1471 | 1471 | 1472 | 4.735 | 4.584 | 6.497 | 53.873 |
| tokyo-run3 | 2405 | 2405 | 2410 | 9.868 | 9.834 | 17.443 | 189.911 |

FLOATは8,258→8,269。検出なしの全6 CSVは2972bf5の全列と完全一致した。
CMCのみも検出なしと既存26列が一致し、両方とDopplerのみも既存列が一致した。
採用解ではDoppler式の超過で6,198衛星アークをリセットし、CMCによるリセットは0。
これは採用解だけの診断集計で、失敗した更新を含む全試行の検出数ではない。
検出数は実スリップの正解ラベルではなく、加速度・雑音との切り分けも残る。

表は各設定の全FLOATを使った、reference.csvをアンテナとみなす参考誤差。
成功時刻が異なり、誤差の変化も一様でないため、精度改善とは判断しない。
IMU仮定とTokyo両レバーアームを含む全仮定・無効理由・検出集計を[report.json](report.json)に保存。
[worst-positions.json](worst-positions.json)には事後調査用の上位10外れと全出力列を残した。

最大外れのNagoya run2週2323/TOW 556760では、4衛星中3衛星だけに履歴があり、
超過0、ambiguity初期化1衛星だった。初回/再捕捉したG32には履歴がなく、
継続観測の検出式ではその時刻の新規衛星の異常を検査できない。
次は初回採用の品質確認と元C++の残差・外れ値処理を比較する。

## 再実行と保存

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_float
python3 -B scripts/run_ppc_slip.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float \
  --output-dir /tmp/new-ppc-slip-comparison

python3 -B docs/comparison/2026-10-21-ppc-slip-detection/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

生CSV・stderr、終了コード、実行コマンドとSHAを保持した。
元C++28ケースは再生成して保存CSVに完全一致。Rust291テスト、examples、fmt、
全ターゲットclippy、releaseが通過し、qualification.jsonと圧縮ログに保存した。
元C++whole RTKは今回再実行せず、先行する既定SPP初期化の記録を参照する。
PPCでのRTK・PPP同等性能は未達。確認済みチェックリストは18/24（75%）。
