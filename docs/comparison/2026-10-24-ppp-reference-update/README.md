# 元C++の精密IF単一更新を直接移植

固定元 `72f3b7c2c4c088dc575286a6f0edf4e407bc499a` の通常の精密IF経路から、
観測行と1回のカルマン更新を移植した。設定値の探索は行っていない。
ライブラリは `PppMeasurementUpdate::ReferenceSinglePass`、PPC exampleは
`-native-update` で明示的に選択する。後者は前段のIF初期化/SPP初期時計も選択する。

元 `ppp_measurement.cpp::formMeasurementEquations` と `ppp_state.cpp::updateFilter` に合わせた点:

- 位置の偏微分は−LOS。補正済み観測の事前位置の幾何で観測行を固定する。
- 最少4衛星、初回はコード観測だけ。継続アークのlock条件を満たしてから位相を使う。
- コード/位相の元prefitゲート、起動時/収束後の閾値、分散floorを使用する。
- 移動時は共分散で標準化した残差の4σ判定。衛星全体の観測行を外し、同じ事前から再判定する。
- 除外しても全状態・ambiguityの次元を保つ。残る衛星が4未満なら更新しない。
- Joseph形式のカルマン更新を1回だけ確定する。別経路のIEKF、5列のコードrank条件、平均NIS、
  コード/位相postfit上限、ZTD非負条件はこの元処理には入れない。

これらは元の別経路の移植であり、既存ゲートの許容差を緩和してテストを通したものではない。
元と同じ20位置の窓・3D偏差0.1 mでprefit用の収束状態を進める。
数値の有限性・SPDの検査はRustで維持する。

## 検証

未修正の元メソッド `formMeasurementEquations` を64回/生成（32条件×gate有無）実行し、
元 `worstStandardizedResidualRow` と元LDLT/Josephの式で比較値を作成した。
2回の生成CSVは全バイト一致。32条件でprefit行、除外衛星、採否、全状態/共分散、NISを照合。
数値は既存の相対1e-10以内。これは元観測組立の実呼出しと更新式の成分比較で、
`updateFilter` / `processEpoch` 自体の新たな呼出しではない。
初回adapterがu8のPRNを文字で出力した失敗を保存し、adapterのみintへ修正した。

最少4衛星、初回コードのみ→次回位相、LLI後のlock再開、衛星除外後も全9 ambiguityを保持する接続テスト、
prefitの厳密境界と元経路が独自NIS/rank/ZTD条件を加えないテストを追加。
Rust456テスト、fmt、clippy、releaseが通過した。旧Nagoya run3 1C/2Lの72列CSVは全バイト一致。
位相未使用のRMSは `Option::None` / CSV空欄とし、0を捏造しない。
新経路の残差診断は固定した事前の線形残差 `r − H(x_post − x_prior)`。
`PppResidualBasis::FrozenPriorLinearizedPosterior` とCSVに明示し、再線形化したpostfitと混同しない。

## PPC結果

全6走行1C/2Wの58,634エポックと、初期化のみでは全点失敗したNagoya run3 1C/2Lを実行。
旧1条件の再演を含む合計69,036行を保存した。
1C/2WのPPP表示は40,277点だが、位相が参加したのは14,666点、コードのみは25,611点。
既存の位相付き12,634点と40,277点を直接比べて継続性の改善とは扱わない。
また旧経路は位相4行以上、今回の「位相あり」は1行以上であり、14,666点との差も同じ採用条件の改善率ではない。

下表の誤差はreferenceがアンテナ位置という仮定での3D誤差(m)。
IMU仮定、Tokyoの両レバーアーム、共通/新規/脱落集合、最悪点は `report.json` に保持。

| 走行・信号 | PPP表示 | 位相あり | コードのみ | 成功集合RMS | 最大 | 旧経路との共通点 | 共通点RMS 旧→今回 |
|---|---:|---:|---:|---:|---:|---:|---:|
| nagoya-run1-1C2W | 6,325 | 5,710 | 615 | 14.799 | 438.186 | 1,603 | 9.951 → 9.912 |
| nagoya-run2-1C2W | 6,490 | 5,674 | 816 | 49.111 | 545.494 | 1,418 | 11.982 → 9.087 |
| nagoya-run3-1C2W | 3,579 | 2,790 | 789 | 31.035 | 177.574 | 681 | 14.072 → 11.473 |
| tokyo-run1-1C2W | 4,738 | 292 | 4,446 | 40.758 | 208.485 | 1,263 | 7.131 → 7.085 |
| tokyo-run2-1C2W | 7,173 | 85 | 7,088 | 18.655 | 198.110 | 1,874 | 4.479 → 6.899 |
| tokyo-run3-1C2W | 11,972 | 115 | 11,857 | 14.753 | 150.969 | 4,773 | 5.156 → 5.830 |
| nagoya-run3-1C2L | 4,433 | 3,640 | 793 | 41.333 | 241.669 | 587 | 22.515 → 23.038 |

500 m超の外れ、負のZTD、Tokyoで位相が継続しない問題を保持する。
解数の増加は精度・実用性能の達成ではなく、全体18/24（75%）、受け入れ未達。
既存のnative全処理CSV2初期化条件×6走行とも比較したが、今回のnative再実行ではない。
位相OSB、SPP種のGPS/Galileo・QC、時計/初期化/履歴が異なり、同じ処理系とは扱わない。
時間はこの実行範囲の記録のみで、速度/メモリの受け入れは行っていない。

## 残る直接移植

予測/Doppler後にGNSS更新が棄却されても、元C++は予測状態と処理時刻、raw観測による
スリップ履歴を進める。Rustはまだ採用時だけ数値時刻を進め、位相履歴もtransactionalである。
新規/未観測ambiguityの保持、SPP時計の毎エポックリセット、初期dt、SPP位置でのZTD初期化、
native測定分散と複数衛星系も含めて、元 `processEpochStandard` の順序を続けて移植する。
独自の閾値調整や新しい測定モデルを先に追加しない。

## 再実行

```sh
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_ppp
python3 scripts/build_ppp_reference_update_fixtures.py --output-dir NEW_NATIVE_DIR
python3 scripts/run_ppp_reference_update.py --output-dir NEW_PPC_DIR
python3 docs/comparison/2026-10-24-ppp-reference-update/verify.py
```

元ライブラリ/Eigen/PPC/精密製品は既存の環境キャッシュを使用する。
将来ソースが進んだ後の保存記録の検証には `verify.py --source-commit RECORDED_COMMIT` を指定する。
fixtureの生成コマンド・出所・元ライブラリhashは `native-generator.json`、
実測のコマンド・全入力/ソースhashは `predeclared-inputs.json` と各 `execution.json` に保存。
