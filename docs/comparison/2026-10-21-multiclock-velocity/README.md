# GPS/Galileoの別時計Doppler速度をPPCで比較

2026-10-21。衛星系ごとの時計ドリフトを推定する速度モデルを選択式で追加した。
全6走行・通常/検出併用についてGPS単独とGPS/Galileoを比較し、24実行を保存。
**速度出力と参考誤差は改善したが、大きな外れが残り、FIX履歴には融合していない。**

## モデル

`VelocityClockModel::PerSystem`は、H=[−LOS,時計_GPS,時計_Galileo,…]を使う。
異なる系の時計ドリフトが同じとは仮定しない。重み付きHouseholder QRで速度と
系ごとの時計range-rateを推定し、全ての時計・交差相関を含む完全な共分散を保持する。
従来4×4フィールドはその先頭主部分行列、単一の時計表示は最初の系を表し、
全時計の名前と値も別に返す。事後自由度は観測数−3−時計数（最少1）。

既定は従来のSingleSystemで、単一系では従来の演算経路を使う。
σ=0.5 m/s、χ²/自由度≤25、未スケール速度分散≤400 (m/s)²、事後スケールの最少1は維持。
ランク不足へ事前拘束や対角下限を加えない。

新raw builderはGPS 1CとGalileo 1C/1Xを明示的に選ぶ。
Galileoでは有効Doppler、有効コード、1C優先の順に選び、他の追尾群も除外記録に保持する。
位相の有無は速度選択に使わない。原値・欠測・追尾コード・LLI、暦toe・送信時刻・回転角は
前回と同じく保持。CSVにも選択追尾・除外理由・系ごとの時計を保存する。
受信機時計の違いをモデル化できても、外れ観測の正しさを保証するものではない。

## 独立検証

[生成器](../../../scripts/upstream_multiclock_velocity.cpp)は、元の測定式・重み・QCを
複数時計へ拡張し、独立したEigenの列pivot QR/LDLTで16ケースを計算する。
**元ライブラリの複数時計速度APIを呼んだ結果ではない。** 単一系は引き続き
[元の実関数11ケース](../2026-10-21-doppler-velocity/README.md)と照合する。

- 2/3時計、2地点、雑音、最少観測・ランク不足・外れ値を照合。
  成否・速度・全時計・完全共分散・χ²/スケール/RMSが相対1e-9以内。
- GPS 3＋Galileo 3では、各系単独は衛星数不足でも、別時計の混合で既知速度を回復。
- 一方の時計だけ変えても、速度・他方の時計・共分散が保たれる。
- 1系の場合は旧経路と完全一致。独立時計の1衛星を加えても速度情報が増えないことを確認。
- 行順序変更、Galileo追尾の選択/欠測・代替群の保持を検証。

最初のsingletonテストに使ったGPS 4方向は、元の分散上限を超える弱い幾何だった。
上限は変更せず、広く分布した8方向の対照で情報不変性を検証した。
6追加を含むRust343テスト、examples、fmt、全ターゲットclippy、releaseが通過。
独立基準は再生成して全バイト一致、コマンド・ソース/ライブラリのハッシュを保存した。

## PPC比較

前回と同じFLOAT位置・同じ入力・同じ閾値で、DopplerだけGalileoと別時計を追加した。
GPS対照の全51列が前回と一致。混合条件も既存の測位・整数診断29列は全行一致し、
FLOATやFIX政策へ速度を戻していない。
分母は1 Hz対応対でFLOATに成功した9,727/9,855点で、元5 Hz全観測の速度出力率ではない。

| 走行 | GPS速度 | 混合速度 | GPS・3D速度誤差RMS m/s | 混合・3D RMS m/s | GPS最大 m/s | 混合最大 m/s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 846 | 856 | 3.688 | 2.247 | 34.504 | 19.889 |
| nagoya-run2 | 1398 | 1425 | 5.799 | 2.235 | 87.310 | 34.811 |
| nagoya-run3 | 908 | 935 | 2.736 | 1.974 | 28.216 | 15.826 |
| tokyo-run1 | 1913 | 2013 | 2.425 | 1.787 | 37.079 | 17.477 |
| tokyo-run2 | 1627 | 1640 | 1.840 | 1.277 | 22.769 | 10.287 |
| tokyo-run3 | 2801 | 2827 | 1.663 | 0.796 | 61.629 | 9.878 |
| 合計 | 9493 | 9696 | — | — | — | — |

表は通常FLOATの全成功速度、referenceをアンテナとみなす参考値。
検出併用では9,616→9,823点。全ての速度失敗・既存無効/未対応エポックも保持する。
通常の共通成功は9,482点、GPSだけ11点、混合だけ214点だった。
共通成功時刻の3D速度誤差RMSも全6走行で低下したが、約35 m/sの最大外れは残る。
異なる成功時刻のRMSだけから優劣を決めない。

両位置仮定とTokyo両レバーアームは[report](report.json)に保持する。
参考ENU速度とレバーアームoffsetの時間微分は[前回と同じ規約](../2026-10-21-rover-velocity/README.md)。
候補選択や閾値へ真値を使わず、採点は実行後だけ。[自由度別の診断](observability-diagnostics.json)も保存。
元C++全RTK/速度ラッパーは今回は再実行していない。
ロック・連続確認・移動履歴、RTK/PPP精度・FIX・速度/メモリの受け入れは未完了。

## 再検証

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_candidates
python3 -B scripts/run_ppc_multiclock_velocity.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_candidates --output-dir /tmp/new-multiclock-velocity
python3 -B docs/comparison/2026-10-21-multiclock-velocity/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

verifyは測位ソルバーを動かさず、ハッシュ・対照・集計と両位置仮定の参考指標を再計算する。
将来ソースが変わった場合、本記録を追加したコミットで検証する。
データ出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)、Copyright (c) 2025 Taro Suzuki。
[MITライセンス](../2026-10-21-gps-galileo-float/input-diagnostics/LICENSE-PPC-Dataset)を保存済み。
次は外れ観測の品質と、速度不確かさを考慮した移動履歴を進める。全体18/24（75%）を維持する。
