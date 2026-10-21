# PPC全6走行の読込と精密製品取得

2026-10-21。PPC実データへの接続を再開した。
以前のCMC実装・比較結果は未復旧だが、入力データのキャッシュは利用できた。
今回の結果は入力準備の検証であり、RTK/PPP解や同等性能の受け入れ結果ではない。

## 入力互換性の修正

変更前の58562c6では、CLIによる全6走行がRINEXの空欄位相シフトで停止した。
[記録](baseline-cli-parse-failures.json)を残した。

- RINEX 3の位相シフト宣言を、空欄・明示ゼロ・非ゼロ、適用衛星と継続行を含めて
  順序付きで保持する。ファイル生成側が適用済みの補正であり、測定値へ再加算しない。
- 観測・航法レコードの空白で桁合わせされたPRN（G 9等）を読めるようにした。
  公開SatelliteIdの通常の文字列書式は変更していない。
- 混合航法データから対象衛星系を明示して選べるAPIを追加した。
  対象外の航法レコードは生データと行番号を保持する。従来の厳密な読込APIも維持する。
- IGSの混合AS/AR時計から衛星ASを選び、受信機ARを値・時刻・名前付きで別保持する
  APIを追加した。受信機時計を衛星時計として使わない。
- Bias-SINEXの自由文コメントを扱い、CASの受信機専用行を未解釈のまま別保持して
  衛星行を選べるAPIを追加した。DSBを絶対OSBへ変換する基準は今回導入していない。

新しいAPIは read_navigation_for_systems、read_satellite_clocks、
read_satellite_biases。GPSだけの入力比較を明示し、他衛星系の実装済みとは扱わない。
CLIの同期処理や測位モデルへの接続はまだ完了していない。

## データと観測時刻の検証

出典は [taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)。
ローカルのZIPと全6走行の観測・航法、およびREADME/LICENSEの内容一致をSHA-256で確認した。
ZIP名からデータのGitコミットを推測せず、ZIPそのものと各入力のハッシュを記録した。
生データの利用条件はPPC-DatasetのMIT License、Copyright (c) 2025 Taro Suzuki。
このディレクトリには生観測を再配布していない。

Rust読込例で全入力を最後まで読み、時刻が厳密に増加することを確認した。
同時刻の候補対を数え、片側だけのエポックも分母から除かず記録した。

| 走行 | 移動局エポック | 基準局エポック | 同時刻の対 | 共通GPS L1 C/Aが4衛星以上の対 |
| --- | ---: | ---: | ---: | ---: |
| nagoya-run1 | 7602 | 1531 | 1521 | 1253 |
| nagoya-run2 | 9451 | 1891 | 1891 | 1411 |
| nagoya-run3 | 5201 | 1041 | 1041 | 702 |
| tokyo-run1 | 11928 | 2391 | 2388 | 1640 |
| tokyo-run2 | 9151 | 1831 | 1831 | 1488 |
| tokyo-run3 | 15301 | 3061 | 3061 | 2430 |

合計で移動局58,634、基準局11,746エポック、同時刻11,733対。
移動局だけ46,901、基準局だけ13エポックがある。
4衛星以上という値はコード・位相の存在と半サイクルを確認した数で、
標高角・航法品質・モデル残差による測位への採用数ではない。
基準局は1 Hz、移動局は5 Hzであり、同時刻だけの処理を5 Hz測位と扱わない。
GPS以外の航法20,049レコードは選択APIで別保持し、その件数を出力した。

## 取得した精密製品

BKG公開ミラーの通常のIGS製品ディレクトリから、PPCの3観測日と前後1日、
計9日分を取得した。日付は2024年の年積算日201〜206、215〜217。

| 種類 | 製品 | ファイル数 |
| --- | --- | ---: |
| 軌道 | IGS0OPSFIN、15分間隔SP3、GPS | 9 |
| 時計 | IGS0OPSFIN、30秒間隔CLK、GPS | 9 |
| コードバイアス | CAS0OPSRAP、日単位DSB | 9 |

全27ファイルでHTTP成功、gzip展開、基本的な書式・日付・データ有無を検証し、
圧縮前後のサイズとSHA-256、出典URLを
[product-acquisition.json](product-acquisition.json)へ保存した。
圧縮合計27,242,490 bytes、展開後113,405,946 bytes。
元の製品は環境内の次のディレクトリへ保存した。

    /workspace/.cloud-setup/datasets/PPC-Precise/igs-20261021

Rustでも9日分すべてのSP3・衛星CLK・衛星DSBを読み込めた。
時計803,519、軌道27,648サンプル、衛星DSB7,542行を確認。
受信機時計513,427レコードと受信機バイアス46,008行は衛星製品と分けて保持した。
各日の開始450秒後だけを補間確認し、各日31衛星で位置・速度・時計を取得した。
残る1衛星/日は利用不能として数え、利用可能扱いにしていない。
全PPC走行時刻の補間・補正の検証は今後行う。

これらは以前要求したCODE MGXのOSB製品ではない。CASのDSBは相対値であり、
IGS時計の信号基準と各観測コードへの適用基準を確認する必要がある。
PPPにはまだ適用していない。受信機アンテナやレバーアームの選択も今回していない。
reference.csvの位置については、以前のユーザー指定どおりIMU中心／GNSS位置の
両仮定で参考評価する方針を維持する。

## 保存・再検証

ソース、入力、バイナリのハッシュを保存してからrelease例を2実行した。
両実行とも終了コード0。生の標準出力/標準エラー、コマンドと終了コードを保存した。
集計は [report.json](report.json)、入力とソースのハッシュは
[predeclared-inputs.json](predeclared-inputs.json)。

追加6テストを含めRust266テスト（docを含む）、examples、fmt、全ターゲットclippy、
releaseビルドが通過した。証拠は [qualification.json](qualification.json) と圧縮ログ。

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_input_audit
python3 -B scripts/run_ppc_input_readiness.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --dataset-zip /workspace/.cloud-setup/downloads/PPC-Dataset-460b21747a78.zip \
  --products-root /workspace/.cloud-setup/datasets/PPC-Precise/igs-20261021 \
  --binary target/release/examples/ppc_input_audit \
  --output-dir /tmp/new-ppc-input-readiness
~~~

製品を再取得するときは acquire_ppc_precise.py に新規の --output-dir を指定する。
保存済み証拠だけの検証は次で行える。ソルバー・ネットワーク・入力キャッシュは不要。

~~~bash
python3 -B docs/comparison/2026-10-21-ppc-input-readiness/verify.py
~~~

次の作業は実運用の時刻同期とキネマティックRTK、PPPの信号補正基準への接続。
両ソルバーのPPC性能受け入れは未達で、全体のチェックリストは18/24（75%）のまま。
