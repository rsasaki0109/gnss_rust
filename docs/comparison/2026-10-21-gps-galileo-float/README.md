# GPS/Galileo FLOATと元C++全RTKのPPC比較

2026-10-21。GPS L1にGalileo E1を加え、衛星系ごとに基準衛星・時計・相関雑音を分けた
FLOATを実装した。全6走行でRust 42実行、元C++全RTK 12実行を保存。
独立デコーダの**345,469入力行が、追尾コード・値・欠測・継続通知を含めて一致**した。

RustのFLOATは8,258→9,727、検出との併用で9,855へ増えた。
精度は走行によって悪化もあり、既定GPSのみを維持する。
**移動/複数衛星系FIXとPPP・速度/メモリの受け入れは未完了。同等性能の受け入れではない。**

## 実装と独立検証

RtkConstellations::GpsGalileoを選ぶとGPS 1CとGalileo E1を使う。
Galileoは各受信機で完全なコード/位相を優先し、1C、次に1Xを選ぶ。原追尾コードは
tracking_pairs()に保持し、選択の変更時はLLIがなくてもその衛星のアークを再初期化する。
GPSは従来どおり厳密に1Cのみ。SPP初期値を選ぶ場合も、そのSPPはGPS 1Cのまま。

DDは衛星系内だけで作り、系間の受信機時計/バイアスを相殺したと仮定しない。
系ごとに参照・Rブロックを持ち、合計3つ以上の独立したDDと位置幾何のランクを要求。
単独1衛星の系は明示的に除外する。SD ambiguityは同じ状態に保持するが、
整数候補の変換も系内だけで行う。既存の初期FIX検証・全基準監査は複数系を明示的に拒否する。

幾何診断はH=[-LOS,時計_GPS,時計_Galileo]の無重みHouseholder QRから計算する。
完全な複数時計共分散と系名も保持し、既存4×4フィールドはその先頭の主部分行列を表す。
事前拘束・対角下限によってランク不足を通す処理はない。

- 1/2/3時計・2地点の24ケースを独立のEigen計算と相対1e-10で照合。
- 元C++物理モデルの移動32エポックで、異なるGPS/Galileo受信機時計と
  Galileo基準局の共通0.25 cycleを与え、全32点が既存3 mm真値許容差を通過。
- GPS 3＋Galileo 3、Galileoだけ4衛星の初期値は、有限の900 m²/900 cycle²事前分布による
  バイアスも含めて元の軌道・DD共分散・Kalman部品から別生成した3初期ケースと比較した。
  位置1 µm、共分散相対1e-10で一致。少数衛星の初期バイアスはそれぞれ約3.865 mm/77.116 mm。
  GPS 3＋Galileo 3は2エポック目以降、3 mm真値許容差を満たした。
- 基準衛星の消失/再捕捉、追尾変更、欠測/半サイクル、系ごとの継続性、
  部分的な観測数不足、重複・スリップ・拒否時の状態不変、未対応FIXの拒否を検証。

部品検証・合成真値・実データの評価を分ける。元RTKProcessor全体の内部同一性は主張しない。

## PPCでの凍結した実行

同じ11,733対応対、仰角10度、基線上限10 km、測量済み基準局、同じ原ファイルを使った。
片側だけの46,901移動局/13基準局エポックも全て出力する。
元5 Hz移動局58,634全時刻の測位率ではなく、以下の率の分母は1 Hzの対応対。
各実行前に入力・ソース・バイナリ・設定を凍結し、真値は事後評価だけに使った。

| 走行 | 対応対 | GPS FLOAT | GPS/Galileo FLOAT | 検出併用FLOAT | 元C++既定SPP RTK | うちFIX |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 1521 | 840 | 861 | 867 | 1374 | 24 |
| nagoya-run2 | 1891 | 1241 | 1438 | 1456 | 1695 | 32 |
| nagoya-run3 | 1041 | 697 | 940 | 956 | 925 | 3 |
| tokyo-run1 | 2388 | 1604 | 2017 | 2048 | 2047 | 27 |
| tokyo-run2 | 1831 | 1471 | 1640 | 1661 | 1675 | 24 |
| tokyo-run3 | 3061 | 2405 | 2831 | 2867 | 2869 | 34 |
| 合計 | 11733 | 8258 | 9727 | 9855 | 10585 | 144 |

FLOAT率は70.4%→82.9%（検出併用84.0%）。元C++既定SPPは90.2%、
FLOAT 10,441 / FIX 144 / SPP退避878 / INVALID 270だった。
SPP退避をRTK成功には数えていない。Rustの複数系パスは全てFLOAT。

元C++の前回RTK位置を優先する別条件も全6走行で保存し、RTKは1,406、
SPP退避10,010、INVALID 317。この設定を実用の比較基準として有利に使わない。
元C++のAR・QCは既定のまま有効。元ソース/静的ライブラリには変更を加えていない。
入力追尾の選択・未対応時刻の継続通知をadapterでそろえても、採用/失敗後の履歴や
SPPの使用衛星/モデルは同一ではない。航法デコード値まで一致したとの主張もない。

| 走行 | Rust GPS/Galileo H RMS m | V RMS m | 最大3D m | 元C++既定SPP H RMS m | V RMS m | 最大3D m |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 2.529 | 7.306 | 142.817 | 6.772 | 19.579 | 373.050 |
| nagoya-run2 | 9.281 | 14.634 | 310.937 | 10.691 | 35.396 | 287.512 |
| nagoya-run3 | 9.147 | 24.874 | 96.814 | 12.612 | 27.140 | 565.993 |
| tokyo-run1 | 10.859 | 19.961 | 141.227 | 12.467 | 23.015 | 537.728 |
| tokyo-run2 | 5.963 | 14.591 | 151.225 | 3.981 | 8.429 | 78.043 |
| tokyo-run3 | 10.029 | 25.910 | 490.732 | 8.306 | 11.243 | 135.073 |

表はreference.csvをアンテナとみなす参考誤差で、それぞれの全RTK出力を使う。
成功時刻が違い、全体の優劣・同等性をRMSの大小だけで判断しない。
IMU位置仮定とTokyo両レバーアームも[root report](report.json)と
[元C++ report](native-comparison/report.json)に保持した。
[共通成功時刻の診断](shared-epoch-diagnostics.json)では、共通/片側/両側不成立の数と
両者の参考誤差も再計算した。元C++のFIX別誤差も保存し、状態名だけで正解と扱わない。
誤FIXの受け入れ閾値や速度/メモリ条件はまだ定めていない。

GPSのみ・GPS＋検出の全CSVは従来の全列と一致した。
追加した5条件はGalileo／検出併用／GPS-SPP初期値／その検出併用／さらに3 mゲート。
SPP＋検出＋ゲートは3,458出力まで減り、これだけを選んで改善とは主張しない。
PPCでの既定設定の変更は今回行っていない。

## 入力診断・ラッパーの確認

[読込専用の診断](input-diagnostics/)では、初期SPPの約1,573 m外れがGPS 4衛星・
残差RMS約2.6e-9 m・GDOP約47で起きたことを確認した。
同じ時刻にはGalileoのコード/位相/基礎航法条件を通る4衛星が利用可能だった。
全6走行で、仰角や位置の検査前の候補はGPS 4衛星以上8,924対、
GPS/Galileoで3 DD以上10,864対。これは入力可能性であり測位成功数ではない。

QZSSは最悪点で位相欠測/非健全航法により補えず、Tokyoの全QZSS選択では
ゼロ軌道レコードを厳密読込が拒否した。この失敗も保存している。
新しい処理はGPS/Galileoを明示的に選び、他系は既存APIの除外記録として保持する。
データを修正したり、航法の健全性・年齢条件を緩めて通したりしていない。

元C++adapterを合成移動32点でも実行し、両初期位置条件とも11 FLOAT / 21 FIX。
最大3D誤差は前回RTK優先約0.865 mm、既定SPP約0.107 mm、入力512行もRustと一致した。
これはラッパーの動作確認であり、PPC性能の代替ではない。

## 保存・再検証

Rust316テスト、examples、fmt、全ターゲットclippy、releaseが通過。
新しい入力台帳helperにも追加のfmt/clippy/release確認を保存した。
元C++の24幾何・32移動・3初期ケースは再生成して全てバイト一致。
生CSV・stderr・両入力台帳・終了コード・コマンド・ハッシュを保持した。
元C++12実行はRust主実験とは別のmanifestで凍結し、元の主実験記録を変更していない。
各manifestのnative実行フラグは、それぞれのrunnerの範囲を表す。

~~~bash
source /workspace/.cloud-setup/activate.sh
cargo build --release --locked --offline -j2 \
  --example ppc_rtk_float --example ppc_mixed_input_ledger
python3 -B scripts/run_ppc_galileo.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float --output-dir /tmp/new-ppc-galileo
python3 -B scripts/run_native_galileo.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --native-build /workspace/.cloud-setup/solver-reference-20261021 \
  --rust-ledger-binary target/release/examples/ppc_mixed_input_ledger \
  --rust-archive /tmp/new-ppc-galileo

python3 -B docs/comparison/2026-10-21-gps-galileo-float/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

証拠のverifyはソルバーを実行せず、両台帳345,469行を再比較する。
将来ソースが変わった場合は、この証拠を追加したコミットで検証する。
データ出典は[taroz/PPC-Dataset](https://github.com/taroz/PPC-Dataset)、
Copyright (c) 2025 Taro Suzuki。[MIT License](input-diagnostics/LICENSE-PPC-Dataset)を同梱。
取得ZIP/入力の出典とハッシュは[入力読込記録](../2026-10-21-ppc-input-readiness/README.md)も参照。

次はSD ambiguity共通オフセットの整合と、移動/複数系FIXを進める。
全対象信号・PPP補正基準・実時間とメモリの比較も未完了。全体18/24（75%）を維持する。
