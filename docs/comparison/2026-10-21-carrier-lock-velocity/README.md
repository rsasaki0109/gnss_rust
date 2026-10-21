# Doppler入力へ位相の有効性条件を加える

2026-10-21。Doppler速度の入力へ、有限の搬送波位相と正常なLLIを要求する選択式条件を追加。
全6走行・通常/検出併用についてGPS/Galileo・条件有無の48実行を保存した。
**RMSは下がるが出力が減り、最大誤差が悪化する走行も残る。既定無効・未融合を維持する。**

## 原観測からの切り分け

[前回の速度比較](../2026-10-21-multiclock-velocity/README.md)の各走行の最悪点を、
読込専用[probe](read-only-probe/report.json)で調べた。
5 Hzで直前の同じ追尾群の位相を保持し、位相差×波長＋現在Doppler×波長×時間差を記録。
系ごとの共通時計を考慮する参考として、選択衛星内の中央値からの差も事後計算した。
全6走行の最悪点で、速度計算に使った一部衛星は現在または直前の位相が欠測していた。
位相を比較できた衛星では、この中央値からの差はおおむねcm単位だった。

probeはソルバーの設定や状態を変更せず、元のsrcは`53f5ded`へ凍結した。
入力ハッシュは前回凍結したファイルと照合した。最悪点の選択は事後の原因調査であり、
独立した未知データでの検証とは扱わない。

## 追加条件

`VelocityEpochConfig::require_carrier_lock`をtrueにすると、選択済みDopplerのrawに
有限のcarrier phaseがあり、LLIのloss-of-lock/half-cycle bitが立っていないことを要求する。
欠測・不正位相・LLI異常は別の除外理由で、原値も保持する。実測0 cyclesは有効。
位相値を速度の式へ入れず、時刻間の整合検査・最小ロック時間・整数の連続確認を
実施したことにはしない。別の追尾群を都合よく選び直す処理も行わない。

既定false、GPS/Galileoの独立時計・σ・χ²/幾何上限は前回のまま。
データ調査を受けた構造的な採用条件であり、数値閾値を真値へ合わせていない。
同じPPCを開発に使った結果であり、hold-out精度検証ではない。

## PPC結果

| 走行 | 混合・無効の速度数 | 有効の速度数 | 無効・3D速度誤差RMS m/s | 有効・3D RMS m/s | 無効最大 m/s | 有効最大 m/s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 856 | 853 | 2.247 | 0.856 | 19.889 | 17.571 |
| nagoya-run2 | 1425 | 1394 | 2.235 | 1.117 | 34.811 | 20.798 |
| nagoya-run3 | 935 | 920 | 1.974 | 1.373 | 15.826 | 15.824 |
| tokyo-run1 | 2013 | 1985 | 1.787 | 0.905 | 17.477 | 25.349 |
| tokyo-run2 | 1640 | 1628 | 1.277 | 0.717 | 10.287 | 11.955 |
| tokyo-run3 | 2827 | 2729 | 0.796 | 0.716 | 9.878 | 16.283 |
| 合計 | 9696 | 9509 | — | — | — | — |

表は通常FLOAT上のGPS/Galileoで、referenceをアンテナと仮定した参考速度誤差。
検出併用では9,823→9,635点、GPSだけでは通常9,493→7,594、検出併用9,616→7,666点になった。
混合系は出力の減少を抑えられるが、位相があるだけでは速度の正しさを保証できない。

通常混合の共通成功時刻でも、3D速度誤差RMSは全6走行で下がった。
共通/片側のみの点数・両位置仮定・Tokyoの両レバーアームは[report](report.json)に保存。
速度を得られなかった点や、元の未対応/無効時刻を隠さず、元5 Hz全観測の出力率とも区別する。
無効24実行の全列は前回と一致し、全48実行で既存の測位/整数診断29列が一致した。
速度はFLOAT・FIXへ戻さず、時刻や整数政策も進めていない。

1追加を含むRust344テスト、examples、fmt、全ターゲットclippy、releaseが通過。
既存の元C++単一時計・独立Eigen複数時計の数値回帰も維持。
元C++全RTK/速度ラッパーの再実行や、RTK/PPPの実用受け入れは今回行っていない。

## 再検証

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_candidates
python3 -B scripts/run_ppc_carrier_lock_velocity.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_candidates --output-dir /tmp/new-carrier-lock-velocity
python3 -B docs/comparison/2026-10-21-carrier-lock-velocity/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

verifyはソルバーを実行せず、ハッシュ・全列/台帳・指標・原観測probeを再確認する。
将来ソースが変わった場合、本記録を追加したコミットで検証する。
データ出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)、Copyright (c) 2025 Taro Suzuki。
[MITライセンス](../2026-10-21-gps-galileo-float/input-diagnostics/LICENSE-PPC-Dataset)を保存済み。
次は残る速度外れと不確かさを調べ、移動体FIXの履歴へ接続する。全体18/24（75%）を維持。
