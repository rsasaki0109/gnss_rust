# GPS/Galileoの整数候補を同じ観測で検証する

2026-10-21。`validate_candidate`をGPS/Galileoの系内DDへ拡張し、
系ごとの全基準衛星の組合せでも同じ整数候補を検証するAPIを追加した。
**読込専用の候補診断であり、移動体FIX政策の実装・FIX出力の増加ではない。**

## 変更

FLOATの不変snapshotが、参加する全ての系の基準衛星を保持する。
整数の順序は、系ごとの基準−対象DDの順序に一致する。
条件付き位置は複数系間の相関を含む完全なjoint covarianceから計算し、
位相/コード残差とcode-minus-phaseはそれぞれの系内だけで評価する。
観測雑音のRも系ごとの参照相関ブロックを保ち、対角化しない。
既存の局所閾値は全て維持した。

`audit_candidate_system_references`は、GPSとGalileoで基準を一つずつ選ぶ
Cartesian productの全組合せへ、同じ整数を系内だけで変換する。
新しいLAMBDA探索や有利な基準の選別は行わない。
callerの上限（1～4096組）を超えれば評価前にエラーとし、不完全な監査を成功と返さない。
整数値と系内の最大spanも、評価前に厳密な2^52未満の範囲を確認する。

従来GPS用の全基準APIは維持する。GPSだけでは新旧の診断結果が完全一致した。
既存`RtkFixPolicy`は、混合系と移動モデルを引き続き拒否し、履歴も進めない。
候補検証を通ることと、ロック・ratio・連続確認・移動履歴を経てFIXを宣言することを区別する。

## 検証

- 元C++の独立物理モデルで作った合成移動32点の既知整数を局所検証し、
  全32点で既存の3 mm真値許容差を満たした。FLOAT状態は不変。
- 4 GPS＋4 Galileoの全16基準組合せで、同じ候補の条件付き位置1 µm、共分散1e-10以内に一致。
- GPSとGalileoそれぞれの誤整数を全16組合せで棄却。
- 組合せ上限・不正次元・系内整数span、公開表示フィールドの改変を検証。
  検証には元の不変snapshotだけを使う。
- 元のGaussian条件付け5ケース、DD基準変換と既存GPS FIX検証も引き続き通過。

6追加を含むRust330テスト、examples、fmt、全ターゲットclippy、releaseが通過。
これは部品・合成の確認であり、元RTKProcessor全体の実データ性能受け入れではない。

## PPCでの読込専用診断

実行前にソース・入力・binary・設定を凍結し、全6走行でGPS/Galileoと検出併用の
12診断を実行した。加えてGPS/Galileo/検出併用の18対照を再実行し、
前回の全CSV列が完全一致した。候補診断でも全時刻のFLOAT位置・使用衛星・reset数が対照と一致。
片側だけの時刻、失敗、候補探索・局所検証のエラーも保存する。
どの入力・候補も真値で選別せず、真値は事後の参考採点だけに使った。

| 条件 | FLOAT snapshot | ratio≥3 | 局所通過 | ratio＋局所通過 | ratio＋全基準通過 | さらに既存4 DD以上 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GPS/Galileo | 9727 | 2347 | 404 | 124 | 61 | 21 |
| 検出併用 | 9855 | 2218 | 517 | 145 | 81 | 25 |

全基準の評価はratioと選択基準の局所検証が通った候補だけで行った。
監査は全て完了し、上限超過で部分監査になった候補はない。
ratio≥3は有限値だけを数えた。最小ロック・連続確認・移動履歴はまだ評価していない。

全基準を通った61/81点のうち40/56点は3 DDだけで、既存の4 DD以上を満たさない。
この診断の局所通過を、そのままFIX数として報告しない。
アンテナ位置仮定で候補の最大3D誤差は、通常GPS/Galileoで約56.531 m、検出併用で約54.049 m。
これらはロック・履歴等を通したFIXではなく、局所検証だけでは不十分であるという結果。
4 DD以上の追加条件も[別の事後診断](eligibility-diagnostics.json)に残したが、FIXは宣言していない。
両位置仮定とTokyoの両レバーアームは[report](report.json)に保持した。

元C++全RTK・入力台帳は今回は再実行していない。
[前回の同入力比較](../2026-10-21-gps-galileo-float/README.md)を参照する。
速度・メモリ、誤FIX、収束・継続率の性能受け入れは未完了。

## 再検証

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 \
  --example ppc_rtk_float --example ppc_rtk_candidates
python3 -B scripts/run_ppc_mixed_candidates.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --float-binary target/release/examples/ppc_rtk_float \
  --candidate-binary target/release/examples/ppc_rtk_candidates \
  --output-dir /tmp/new-ppc-candidates
python3 -B docs/comparison/2026-10-21-mixed-candidate-validation/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

verifyはソルバーを再実行せず、ハッシュ・対照との一致・集計を再計算する。
将来ソースが変わった場合、この記録を追加したコミットで実行する。
データ出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)、Copyright (c) 2025 Taro Suzuki。
[保存済みMITライセンス](../2026-10-21-gps-galileo-float/input-diagnostics/LICENSE-PPC-Dataset)を参照。

次はロック・整数の連続確認と、移動体に適した履歴の接続を進める。
部分AR・FIX維持、PPP、実用性能の受け入れも残り、全体は18/24（75%）を維持する。
