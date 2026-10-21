# 全基準衛星で既存の局所ゲートを要求するRTK FIX方式

2026-10-21。[直前の診断](../2026-10-21-fix-reference-audit/README.md)で、
同じ候補の条件付き位置・共分散が一致しても、位相RMSによる判定が基準衛星に
依存することを確認した。今回、その候補が全基準衛星の局所検証を通った場合に
連続確認へ進める、選択式の保守的なFIX方式を追加した。

## 実装

APIは RtkFixPolicy::with_reference_gates(config, ReferenceGateMode::AllReferences)、
CLIは rtk の --all-reference-gates。従来のコンストラクターとCLIの既定動作は継続する。

元基準のLAMBDA候補を一度だけ選び、整数差を厳密に変換して全基準で再検証する。
全結果の共通通過を要求し、棄却理由は既存のゲート順で最初のものを返す。
最良基準の選択、閾値の変更、条件付き値のFLOATへのフィードバックは行わない。
いずれかの数値評価が失敗した場合も採用しない。

不採用時は連続確認を解除し、既存のロック・ratio・再確認・静止履歴を適用する。
reference_auditに全基準の判定を保持し、validationは元基準での検証結果を保持する。
そのため全体の採否は RtkFixAttempt.decision / fixed_solution で確認する。

これは同じ候補と観測に対する局所ゲートの対策であり、
基準衛星ごとのLAMBDA候補選択・ratioまで同じになるという主張ではない。
全基準の条件付き位置と残差を評価する追加計算があり、処理性能は未評価。

## 比較条件と保存物

既存のC++生成済みGPS L1合成RINEXを用い、次の4条件を固定した。

- clean: 雑音なし、32エポック・9衛星。
- noisy: コード・位相に元の固定雑音がある同数の入力。
- phase-bias: cleanのG01 L1Cに0.15 cycle、ゼロ始まり10〜31番目で加算。
- transient: cleanのG01 L1Cに0.25 cycle、6番目だけ加算。

生成した2ファイルを保存し、数値フィールド以外（時刻・LLI等）は変更していない。
固定幅の末尾フラグ欄を維持するため、Gitの末尾空白検査はこの2入力だけ除外している。
全Rustソース、入力、生成器とテスト、バイナリのSHA-256を
[predeclared-inputs.json](predeclared-inputs.json)へ保存してから、
release版CLIで各条件を従来方式／全基準方式の双方で実行した。
以後のソース・入力変更はない。

各実行のコマンド・終了コード・非圧縮SHA-256と、生のCSV/標準エラーのgzipを保存した。
[report.json](report.json)は全32エポックを時刻で照合し、欠測やFLOATを除外せず集計した。
新しいC++全ソルバー比較、PPC実データの測位、実時間・メモリの受け入れ計測はない。

## 実測した合成ケースの結果

8実行すべて終了コード0、計256エポック行を確認した。

| 条件 | 従来方式のFIX /32 | 全基準方式のFIX /32 | 全基準方式で位相残差棄却 |
| --- | ---: | ---: | ---: |
| clean | 27 | 27 | 0 |
| noisy | 27 | 27 | 0 |
| phase-bias | 24 | 5 | 22 |
| transient | 25 | 25 | 1 |

phase-biasでは偏りが始まる前の5 FIXだけを維持し、偏りのある22エポックを
すべてFLOATへ戻した。整数真値が正しくてもモデル不整合なら棄却する保守方式で、
このFIX減少を実データでの誤FIX率改善とみなす評価は行っていない。
transientでは6番目を棄却し、7番目で再確認、8番目にFIXへ復帰した。

全条件・全エポックでFLOAT更新NIS、参照衛星、衛星数、反復数、ambiguity reset、
候補ratioと組数は両方式で一致した。共通FIXエポックの出力座標差は0 m
（CLIの小数6桁出力）。状態不変性は下記の単体/統合テストで別途検証した。
最大残差列は丸めた表示値であり、実際のゲートは丸め前の値を使う。

## 回帰検証

追加6テストで以下を確認した。

- 正常2条件で従来の27 FIXと位置・候補整数が一致し、FLOAT状態を変更しない。
- 偏りのある22エポックの棄却、一時異常後の再確認、政策reset後のモード維持。
- 実際の不変スナップショットの元基準を9通りに変更しても、同じ候補の
  全基準結果と集約判定が一致する。正常／偏りと、元候補／1整数を誤らせた候補で確認。
- CLIの追加統計、未評価時の空欄、FIX専用オプションと重複指定の検査。

[qualification.json](qualification.json)と圧縮ログにRust260テスト（docを含む）、
examplesのコンパイル、fmt、全ターゲットclippy、releaseビルドの成功を保存した。
以前の254テストに加えて6テスト。既存の誤整数テストにも集約棄却の確認を追加した。

## 再実行

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2
python3 -B scripts/compare_rtk_reference_gates.py \
  --binary target/release/gnss-rust \
  --output-dir /tmp/new-all-reference-gates
~~~

出力先には新規ディレクトリが必要。
保存済み結果と当時のソース・入力の整合性だけを調べるには次を使う。
GNSSソルバーやネットワークは実行しない。

~~~bash
python3 -B docs/comparison/2026-10-21-all-reference-gates/verify.py
~~~

過去の証拠ディレクトリは変更していない。過去の検証器は当時のソースのSHA-256を
要求するため、それぞれの保存コミットで実行する。
PPCでのRTK・PPP実用同等性能は未達で、確認済みチェックリストは18/24（75%）。
