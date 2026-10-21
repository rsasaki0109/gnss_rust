# PPC移動局のGPS Doppler速度を確認する

2026-10-21。既存GPS/Galileo FLOATの各snapshot位置で、GPS 1CのDoppler速度を計算した。
全6走行で通常/検出併用の12実行を保存。**FIXや移動履歴への融合は行っていない。**

前回のFLOAT位置・時刻・整数候補・局所検証・全基準監査を含む29列が全行一致した。
速度は追加の診断であり、既存結果へ影響しない。
仰角10度、σ=0.5 m/s、最少4衛星、χ²/自由度≤25、未スケール速度分散≤400 (m/s)²は
[前回の基盤](../2026-10-21-doppler-velocity/README.md)のまま。閾値や信号選択を真値で調整しない。

## 対象時刻と参考速度

1 Hzの11,733対応対のうち、FLOATに成功した9,727/9,855点で速度を試みる。
未対応5 HzエポックやFLOAT失敗時刻に、補間した位置・速度を作らない。
全ての未出力・失敗・速度エラーも元の時刻台帳に残す。

reference.csvに与えられているEast/North/Up VelocityをECEFへ回し、事後の参考採点をする。
位置の基準点が不明という既存の回答に従い、アンテナ仮定とIMU仮定を両方保持する。
IMU仮定には、既存のroll/pitch/headingから回したレバーアームoffsetの時間微分を加える。
通常は前後0.2 sによる中央差分、端点は片側差分、差分区間が1 s超のときは補わない。
Tokyoの両レバーアームも別々に評価する。これは参考速度であり新たな真値基準点の確定ではない。

## 結果

| 走行 | 通常FLOAT上の速度/試行 | 検出併用の速度/試行 | 通常・水平速度誤差RMS m/s | 鉛直RMS m/s | 最大3D速度誤差 m/s |
| --- | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 846 / 861 | 852 / 867 | 1.689 | 3.279 | 34.504 |
| nagoya-run2 | 1398 / 1438 | 1415 / 1456 | 1.808 | 5.510 | 87.310 |
| nagoya-run3 | 908 / 940 | 924 / 956 | 1.205 | 2.456 | 28.216 |
| tokyo-run1 | 1913 / 2017 | 1940 / 2048 | 1.339 | 2.021 | 37.079 |
| tokyo-run2 | 1627 / 1640 | 1648 / 1661 | 1.174 | 1.417 | 22.769 |
| tokyo-run3 | 2801 / 2831 | 2837 / 2867 | 0.822 | 1.446 | 61.629 |
| 合計 | 9493 / 9727 | 9616 / 9855 | — | — | — |

誤差の表はreferenceをアンテナと仮定した参考値。両位置仮定の全指標は[report](report.json)。
通常234、検出併用239点の速度失敗を残した。衛星数不足、残差超過、弱い幾何を区別する。
速度出力率は元の全5 Hz観測58,634点の出力率ではない。

報告された最大軸の標準偏差の中央値は約1.23～1.93 m/sで、基地局より不確かだった。
Nagoya run2の最大約87.31 m/s外れはTOW557223、GPS Doppler 5衛星で起きた。
[衛星数別の事後診断](observability-diagnostics.json)では4衛星だけでなく、5衛星以上にも
大きい外れが残った。単純に4衛星を除くだけで解消したとは扱わない。
推定されたσや局所残差だけを根拠に、この速度を信頼したFIX履歴へ入れない。
次にGalileoの追加と系ごとの時計ドリフト、観測品質を検証する。

Rust337テスト、examples、fmt、全ターゲットclippy、releaseが通過。
速度推定の基盤は前回の元C++11ケースを参照する。今回は元C++全RTK/速度ラッパーを再実行していない。
RTK・PPPの同等性能、FIX率/誤FIX、速度/メモリの受け入れも未完了。

## 再検証

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_candidates
python3 -B scripts/run_ppc_rover_velocity.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_candidates --output-dir /tmp/new-rover-velocity
python3 -B docs/comparison/2026-10-21-rover-velocity/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

verifyはソルバーを再実行せず、入力/ソース/CSVのハッシュ・既存列の一致と指標を再計算する。
将来ソースが変わった場合、この記録を追加したコミットで検証する。
データ出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)、Copyright (c) 2025 Taro Suzuki。
[保存済みMITライセンス](../2026-10-21-gps-galileo-float/input-diagnostics/LICENSE-PPC-Dataset)を参照。
全体は18/24（75%）を維持し、未融合の速度出力をFIX達成に数えない。
