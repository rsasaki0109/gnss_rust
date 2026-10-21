# SD ambiguityの基準整合とPPCでの切り分け

2026-10-21。元C++の系ごとのSD ambiguity共通オフセット処理を選択式で移植した。
独立した元の`updateBias`実呼出しと照合し、PPC全6走行で7条件の有効・無効、計84実行を保存。
**大きな外れ値の主因ではなかった。既定無効を維持し、同等性能の達成とは扱わない。**

## 処理と適用範囲

各衛星の初期値は、SD carrier cycles − SD code metres / wavelength。
継続中の衛星だけから、その初期値と前回状態との差を系ごとに平均し、継続状態へ同じ量を加える。
新規/再捕捉・既知スリップ・追尾変更した衛星は平均に含めず、現在の初期値を使う。
これにより古い状態同士のDDを保ちつつ、新しく加わる状態の基準をそろえる。

`RtkFloatConfig::align_ambiguity_datums`をtrueにすると有効。既定はfalse。
GPS L1/Galileo E1を対象とし、異なる系の平均は混ぜない。
全体リセット・失敗後・長い間隔では継続状態がないため、平均はゼロ。
共分散・雑音・事前分散・スリップ検出閾値・QC閾値は変更しない。
今回の処理は時計推定や追加の独立測定ではない。
候補状態だけで処理し、更新に失敗した場合は数値状態を確定しない。

元実装は数値0を未初期化の印として使うが、Rustでは`Option`で区別する。
**継続中の0 cyclesを捨てない**点は意図した差であり、元のゼロケースも別に記録した。
未対応FIX政策や移動体FIXを有効にする変更は含まない。

## 元の実処理と合成検証

[生成器](../../../scripts/upstream_rtk_datum.cpp)は、固定コミット
`72f3b7c2c4c088dc575286a6f0edf4e407bc499a`の既存friendを使い、未変更の
`RTKProcessor::updateBias`を呼ぶ。元ソース・静的ライブラリは変更していない。
処理だけを分離する部品ケースでは検出器とprocess noiseを無効にする。
PPC実験の検出器や閾値を無効化して同等性を主張するための設定ではない。

- 新規のみ、継続のみ、新規の参加、継続1衛星、2系、1衛星スリップ、全スリップ、
  平均0、大きい共通オフセットの9ケースで元の値と相対1e-12以内に一致。
- 保持した非対角共分散、初期化した900 cycle²と独立性も元の結果で確認。
- ゼロの印の挙動は10番目のケースとして記録。Rustの有効なゼロとの違いをテストする。
- 合成移動32点では、GPSとGalileoに別々の変動する共通位相と消失/再捕捉を与えた。
  無偏りの対照と位置1 µm・共分散相対1e-10以内で一致。
- ゼロ基線/ゼロ状態、追尾変更の平均からの除外、棄却時の状態維持と
  同時刻の正常入力で再試行できることを確認。

最初の失敗回復テストでは、移動2エポック後の入力が元の3 m prefit gateを超えた。
閾値は緩和せず、時刻を進めない失敗の後に同時刻の正常観測を再試行するテストへ修正した。
新規8テストを含む全324テスト、examples、fmt、全ターゲットclippy、releaseが通過した。

## 変更前の読込専用診断

[read-only-probe](read-only-probe/)は変更前のソース（`906f121`）でGPS/Galileo各6走行を実行。
全エポックの状態・位置は前回のCSVと一致し、基準値は読まず、継続中の旧状態と
現在の生コード/位相による初期値だけから差を計算した。
新規/再初期化と継続が混在した行はGPS 284、GPS/Galileo 792系行。
平均差が距離換算1 mを超えたのはそれぞれ25、101系行、10 mを超えたのは後者の2行だった。
成功時だけの診断であり、失敗した更新の差は推定していない。

Nagoya run2、GPS週2323/TOW556760の旧最大外れでは、継続3衛星＋新規G32、
平均差は−0.3133 mだった。約938 mの誤差をこの差だけで説明する証拠はない。

## PPC全6走行の結果

前回と同じ入力ハッシュ、11,733対応対、仰角10度・基線10 km、測量済み基準局を使った。
設定・ソース・バイナリを実行前に凍結し、参考真値は事後の採点にだけ使った。
GPS / GPS＋検出 / Galileo / Galileo＋検出 / Galileo＋GPS-SPP /
その検出併用 / さらに3 m prefit gateの7条件を、有効・無効で実行した。
片側のみや失敗時刻も全て残した。

| 条件 | 無効FLOAT | 有効FLOAT |
| --- | ---: | ---: |
| GPS | 8258 | 8258 |
| GPS＋検出 | 8269 | 8269 |
| GPS/Galileo | 9727 | 9727 |
| GPS/Galileo＋検出 | 9855 | 9855 |
| GPS/Galileo＋GPS-SPP | 9705 | 9705 |
| GPS/Galileo＋GPS-SPP＋検出 | 9834 | 9834 |
| 上記＋3 m prefit gate | 3458 | 3453 |

無効の42実行は、前回の全CSV列と一致した。
通常GPS/Galileoの水平RMS変化は全6走行で最大約1.3 mmだった。
検出併用などでは悪化もあり、精度改善や既定の変更を正当化する結果ではない。
同じ出力時刻での位置差は、通常GPS/Galileoで最大約0.379 m、検出併用で約3.934 m。
集計RMSの変化が小さいことを、全エポックの位置が同一との意味には使わない。
GPSの旧最大3D誤差は937.912→937.882 mで、大きな外れが残った。
ゲート併用の出力は計5点減り、一部のRMS低下だけを改善と扱わない。

両位置仮定とTokyoの両レバーアームは[report](ppc-comparison/report.json)に保持した。
[共通時刻の診断](shared-epoch-diagnostics.json)では、共通/片側/両側不成立の数と
対応する参考誤差、2つの旧外れ時刻も記録した。
元C++の既定SPP出力との比較は[前回の不変記録](../2026-10-21-gps-galileo-float/README.md)を
再利用して明示的に区別した。今回、全native RTKの再実行や入力台帳の再生成はしていない。
速度・メモリや誤FIXの受け入れも未実施。

## 再検証

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 --example ppc_rtk_float
python3 -B scripts/run_ppc_datum_alignment.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float --output-dir /tmp/new-ppc-datum
python3 -B docs/comparison/2026-10-21-ambiguity-gauge/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

verifyは測位ソルバーを再実行せず、入力・ソース・生CSV/ログのハッシュと指標を再計算する。
変更前診断のsrcハッシュは記録したGitコミットから検証する。
将来ソースが変わった場合、本記録を追加したコミットで検証する。
元C++の部品生成コマンドとライブラリ/ソースのハッシュは[component-generation.json](component-generation.json)。
データ出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)、Copyright (c) 2025 Taro Suzuki。
[MITライセンス](../2026-10-21-gps-galileo-float/input-diagnostics/LICENSE-PPC-Dataset)を保存済み。

次は移動体/複数系FIXと初回観測の品質処理を進める。PPP補正基準への接続も残る。
全体は18/24（75%）を維持し、部品移植を実用性能の受け入れに数えない。
