# 元C++の精密IF初期化を直接移植

`ppp_state.cpp::initializeAmbiguityState` に合わせ、新規・スリップ後のIF ambiguityを
補正済み位相 − 予測距離、分散25 m²で初期化する経路を追加した。
従来のRustは独立な平均0・分散1e14 m²で、元実装とは異なっていた。
`initializeFilter` の順序も合わせ、明示した同時刻コードSPPの時計を初期平均へ取り込んでから位相を初期化する。
SPP時計の毎エポック融合・時計ドリフトの取込は行わない。

固定元コミット `72f3b7c2c4c088dc575286a6f0edf4e407bc499a` の実メソッドを32ケースずつ2回実行。
fixtureの全バイトが一致し、Rustも全32ケースの平均・分散が一致した。
生成元、ライブラリ、ヘッダーのSHA-256と実コマンドは `native-generator.json` に保存。
元ソースの修正はない。これは初期化メソッドの比較であり、全PPPProcessorの同等性を示さない。

新規9アーク、健全な継続で再初期化しないこと、LLIによる1アークのリセット、
更新棄却時の提案状態と不正時計の扱いを接続テストで確認した。
Rust451テスト、fmt、clippy、releaseが通過。ログと `checks.json` を同梱。

PPC Nagoya run3、1C/2Lの全5,201エポックを `-native-if-prior` で実行した。
生CSV・stderr・コマンド・ソースハッシュ・状態数は同梱ファイルに保存。
全5,201点が棄却された（PostfitRejected 2,264点など）。この失敗を保持する。
元の初期化とRust独自の更新/ゲートを混ぜた状態で、精度改善や性能同等の受け入れではない。
旧IEKF・測定ゲート・棄却後の数値状態/位相履歴の処理はまだRust独自経路のまま。
次は元の単一カルマン更新とエポック処理順を直接移植する。

ライブラリは `PppAmbiguityInitialization::ReferencePrecise` と
`PppMotionSeedPolicy::initialize_clock_from_seed` で明示的に選択する。
PPC exampleでは `-seed-position` または `-seed-doppler` と `-native-if-prior` を指定する。
全体のチェックリストは18/24（75%）、RTK/PPPの実用同等性能は未達。
