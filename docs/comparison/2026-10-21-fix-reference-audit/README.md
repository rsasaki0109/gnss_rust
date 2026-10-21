# 同じ整数候補の基準衛星別FIX局所判定

2026-10-21。前回の[DD基底の部品診断](../2026-10-21-reference-basis/README.md)を、
実際のRTK FLOATが保持する観測・単差状態・共分散と既存の局所判定へ接続した。
新しい audit_candidate_references APIは各基準衛星での条件付き位置、共分散、
残差統計と最初の棄却理由を返す。元のFLOAT状態とFIX政策の履歴は変更しない。

## 入力と実行

既存のC++生成済みGPS L1合成RINEXを使用した。元C++全ソルバーの新実行、
PPC実データの実行、速度・メモリの受け入れ評価は今回行っていない。
入力の生成元と真値は[fixtureの説明](../../../tests/fixtures/README.md)を参照。

32エポック・9衛星のclean/noisyと、cleanへG01の位相を0.15 cycle加えた
phase-biasの3条件を固定した。偏りはゼロ始まりの10番目以降で、
LLIなどの入力フラグは変えない。条件・入力・主要ソース・実行バイナリのSHA-256を
[predeclared-inputs.json](predeclared-inputs.json)へ保存してから3実行を開始した。
保存後のソース・入力変更はない。

各エポックで元の基準衛星のLAMBDA第1候補を一度だけ求め、全基準衛星へ
同じ候補を整数差分で変換する。基準衛星ごとの再探索・最良基準の選択はない。
合成整数真値は候補の正誤を報告するためだけに用いる。
各入力整数と全基準間の整数差は絶対値2^52未満に制限し、
一つでも表現できなければAPI全体が失敗する。

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example rtk_candidate_reference_audit
python3 -B scripts/run_rtk_candidate_reference_audit.py \
  --binary target/release/examples/rtk_candidate_reference_audit \
  --output-dir /tmp/new-rtk-candidate-reference-audit
~~~

出力先は新規ディレクトリが必要。各実行のコマンド・終了コードと、
変更していない標準出力CSV/標準エラーのgzip、非圧縮SHA-256、
全エポックの基準衛星別判定を保存している。

## 観測結果

全3実行が終了コード0。各32エポック・288基準行を照合した。

| 条件 | 基準衛星で局所判定が異なるエポック | 元基準の局所通過 | 真値と一致した候補 |
| --- | ---: | ---: | ---: |
| clean | 0/32 | 32/32 | 32/32 |
| noisy | 0/32 | 32/32 | 32/32 |
| phase-bias | 22/32 | 30/32 | 32/32 |

phase-biasの差は全て既定の位相残差ゲートによる通過／棄却の差だった。
最初の差はGPS週2300、週内秒345900で、G01基準の位相RMSは約0.02040 m、
G02基準では約0.007918 m。既定のRMS上限0.01 mをまたぐ。
同じ物理候補の条件付き位置の基準間差は全3条件を通じて最大約1.03e-14 m、
条件付き共分散の成分差は約9.64e-16 m²だった。
位置・共分散と相関を保った重み付きスコアが一致しても、単純RMSは一致しない。

ここでの通過数は局所判定だけであり、FIX出力数ではない。ロック・ratio・連続確認・
静止履歴の政策判定をこの実行例では行っていない。
postfit_nis_per_rowという既存APIの列は公称Rによる事後残差の二次形式/行数であり、
校正されたカイ二乗確率や新たなFIX保証として扱わない。

## 検証と次の判断

追加した4テストでは、雑音なし/あり全576基準の独立SD整数真値との照合、
正定値共分散・条件付き位置、誤整数、公開表示フィールドの改変、
次元/厳密整数範囲の境界、固定位相偏りを検証した。
診断あり/なしでFLOAT解・FIX政策の結果と履歴が一致することも確認した。

[qualification.json](qualification.json)と圧縮ログに、Rust全254テスト
（docを含む）、全examplesのコンパイル、fmt、全ターゲットclippy、
releaseビルドの成功を保存した。examples自体には単体テストはなく、
上記の実バイナリ3実行を別途保存している。

現在のソースと保存したログ・CSV・報告の整合性だけを再検証するには次を使う。
GNSSソルバーやネットワーク処理は実行しない。

~~~bash
python3 -B docs/comparison/2026-10-21-fix-reference-audit/verify.py
~~~

既存の閾値やFIX政策は変更していない。次は基準衛星に依存する棄却の扱いを、
誤FIXの抑制と両立する条件で評価する必要がある。
PPCでのRTK・PPP実用同等性能は未達。現在のチェックリストは18/24（75%）のまま。
