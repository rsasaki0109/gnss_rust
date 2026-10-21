# 単一衛星系のDoppler速度とPPC固定基地局での診断

2026-10-21。移動体FIXの履歴に使う速度の基盤として、元C++の
`spp_velocity::solveVelocity`に対応する単一衛星系のDoppler速度・時計ドリフト推定を追加した。
**RTKへの融合・移動体FIX政策への接続はまだ行っていない。**

## モデルと品質検査

準備済み観測のH=[−LOS,1]、測定range rate=−Doppler×c/fを使い、
衛星速度・時計ドリフトを差し引いて、受信機速度3成分と時計range-rateを推定する。
標高角依存の重み付きHouseholder QRを使い、位置の事前分布は使わない。
完全な速度/時計共分散と相関、各観測残差、χ²/自由度を返す。

元の既定と同じσ=0.5 m/s、最少4衛星、χ²/自由度≤25、未スケールの速度分散≤400 (m/s)²。
自由度はmax(1,n−4)、共分散の事後スケールはmax(1,χ²/自由度)とし、
少数衛星のほぼゼロ残差で過信しない。ランク不足に事前拘束や対角下限は加えない。
単一時計の初期実装なので、異なる衛星系の混合は明示的に拒否する。

元はDoppler数値0を欠測の印として除外するが、Rustの準備済み観測では実測0を有効に扱う。
欠測はrawの`Option`で区別する。不正値・重複・設定も明示的にエラーとし、
元の設定値への暗黙の置換や不正観測の黙った除去はしない。

GPS 1Cのraw builderは、2回の送信時刻/時計評価と受信ECEFへの位置・速度の両方の回転を行う。
元と同じ回転規約を使い、独自のΩ×r項を追加しない。
健全性・絶対時刻での暦年齢・仰角を検査する。コード欠測時は元の0.075 s初期値を使い、
実測0/不正コードとは区別する。原追尾コード・値・欠測・LLI等はraw cloneに保持し、
送信時刻・使用暦toe・回転角と、全除外理由も記録する。
位相や整数値は速度推定に使わず、非ゼロの未適用epoch時計オフセットは拒否する。

## 独立検証

固定した元コミット`72f3b7c2c4c088dc575286a6f0edf4e407bc499a`の未変更の
`solveVelocity`を[生成器](../../../scripts/upstream_velocity.cpp)から呼び出した。
静止/移動、雑音4条件、4衛星、観測数不足、弱い幾何、外れ値、ランク不足の11ケースで、
成否・速度・時計・共分散・RMSを照合した。速度/共分散の相対許容差1e-9、
時計ドリフト絶対1e-16 s/s、時計分散1e-23以内。
12番目は元のゼロ欠測印を記録し、Rustが有効な実測ゼロを保持する意図した差をテストする。
元の生成器を再実行してfixtureの全バイト一致を確認し、コマンド・静的ライブラリのハッシュを保存。

7追加テストで閾値ちょうど/1 ULP下、実測ゼロ、コード/位相欠測、時刻・古い暦、
追尾コード・LLIの保持、重複・独立時計の拒否を確認した。
Rust337テスト、examples、fmt、全ターゲットclippy、releaseが通過した。
raw builder全体の元ラッパーとの同一性や、移動局での速度精度はまだ主張しない。

## PPC基地局の全エポック

同じハッシュのbase.obs/base.nav、測量済み基地局位置を使い、GPS 1C・仰角10度・
最大暦年齢14,400 s、上記既定QCで実行した。移動局やreference.csvは読まない。
全時刻キーは前回のPPC入力台帳と一致した。

| 走行 | 全時刻/速度出力 | 静止仮定の速度RMS m/s | P95 m/s | 最大 m/s |
| --- | ---: | ---: | ---: | ---: |
| nagoya-run1 | 1531 / 1531 | 0.01814 | 0.03314 | 0.05253 |
| nagoya-run2 | 1891 / 1891 | 0.01843 | 0.03439 | 0.06213 |
| nagoya-run3 | 1041 / 1041 | 0.01782 | 0.03134 | 0.09662 |
| tokyo-run1 | 2391 / 2391 | 0.05786 | 0.10363 | 0.17843 |
| tokyo-run2 | 1831 / 1831 | 0.05470 | 0.09480 | 0.18252 |
| tokyo-run3 | 3061 / 3061 | 0.07076 | 0.12785 | 0.25603 |

合計11,746点を全て保存した。固定基地局の真の速度をゼロとみなす診断であり、
移動局の5 Hz速度、変位積分、FIXの誤り/継続、RTK・PPPの同等性能や処理性能の受け入れではない。
閾値を調整してこの結果を選んだり、速度をFLOATや整数状態へ入れたりしていない。
既存RTK/PPPコードの変更はモジュール登録だけで、既存の回帰テストも通過した。

## 再検証

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_base_velocity
python3 -B scripts/run_ppc_base_velocity.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_base_velocity --output-dir /tmp/new-ppc-base-velocity
python3 -B docs/comparison/2026-10-21-doppler-velocity/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

verifyはソルバーを実行せず、ハッシュ・台帳・保存済みCSVの指標を再計算する。
将来ソースが変わったら本記録を追加したコミットで検証する。
データ出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)、Copyright (c) 2025 Taro Suzuki。
[MITライセンス](../2026-10-21-gps-galileo-float/input-diagnostics/LICENSE-PPC-Dataset)を保存済み。
次は移動局のDoppler評価と移動履歴への接続を進める。全体18/24（75%）を維持する。
