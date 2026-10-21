# 元C++のSPP初期化を既定へ戻したPPC比較

2026-10-21。[前回の比較](../2026-10-21-native-rtk-common-inputs/README.md)から
元C++のprefer_rover_position_seedだけをtrue→false（既定値）に戻した。
全6走行のGPS 1C・11,733対応対を再実行し、初期位置の選び方による影響を切り分けた。
PPCでの実用RTK・PPP同等性能の受け入れではない。

元ライブラリ・AR・QC・仰角・基線上限・時刻対応付け・LLI通知は変更していない。
wrapperが与える最初のヘッダー位置と以後の前回RTK位置も同じ方針だが、
既定方式ではSPPが失敗した場合の代替としてのみ使われる。
真値やreference.csvは事後の参考評価にのみ使った。

生成した[adapter](native_ppc_rtk_default_seed.cpp)は、元adapterに対して
上記の一つの設定値とstderrの設定説明だけを変更したもの。
元ライブラリのソースを変更したものではない。
生成・構築・実行は[scripts/run_native_ppc_seed_comparison.py](../../../scripts/run_native_ppc_seed_comparison.py)
で行い、実行前に入力・ソース・ライブラリ・バイナリのSHAを凍結した。

## 結果

| 走行 | 前回RTKを優先したRTK出力 | 既定SPPのRTK出力 | うちFIX | SPP退避 | INVALID |
| --- | ---: | ---: | ---: | ---: | ---: |
| nagoya-run1 | 353 | 1124 | 2 | 51 | 346 |
| nagoya-run2 | 222 | 1055 | 3 | 151 | 685 |
| nagoya-run3 | 86 | 509 | 5 | 147 | 385 |
| tokyo-run1 | 273 | 1196 | 1 | 353 | 839 |
| tokyo-run2 | 181 | 1152 | 22 | 240 | 439 |
| tokyo-run3 | 244 | 2002 | 1 | 479 | 580 |
| 合計 | 1359 | 7038 | 34 | 1421 | 3274 |

全11,733対中RTK出力は11.6%→60.0%、SPP退避は5,765→1,421へ変わった。
新しい内訳はFLOAT 7,004 / FIX 34。元5 Hz観測を全て測位した継続率ではない。
同じwrapperで初期位置設定だけを変えた結果なので、この設定変更の影響は確認できる。
Rustとの全機能の同等性や、SPPを常に信用できることを示すものではない。

| 走行 | 前回水平RMS m | 既定SPP水平RMS m | 既定SPP鉛直RMS m | 既定SPP最大3D m |
| --- | ---: | ---: | ---: | ---: |
| nagoya-run1 | 74.363 | 1.358 | 3.829 | 102.684 |
| nagoya-run2 | 79.781 | 14.636 | 69.194 | 462.605 |
| nagoya-run3 | 181.405 | 6.514 | 26.433 | 172.720 |
| tokyo-run1 | 1180.363 | 17.340 | 35.243 | 767.587 |
| tokyo-run2 | 47.694 | 2.519 | 4.766 | 78.922 |
| tokyo-run3 | 46.371 | 4.121 | 4.378 | 70.769 |

表は各設定の全RTK出力を使った、reference.csvをアンテナとみなす参考誤差。
成功時刻が違うので、RMSの比をそのまま精度改善率とはしない。
IMU仮定とTokyo両レバーアームも[report.json](report.json)に保存した。
[diagnostics.json](diagnostics.json)に三者共通のRTK成功時刻だけの参考誤差と
今回のFLOAT/FIX別誤差も保存。今回の34 FIXにも最大約2.93 mの3D誤差があり、
FIX状態だけで正解とは判断できない。誤FIX受け入れの閾値はまだ定めていない。

入力観測の台帳SHAは全6走行で前回に完全一致した。
値と欠測フラグを含む両デコーダの原台帳は前回のディレクトリに保持されている。
全エポックの時刻・片側入力・両末尾は同じで、有効出力の時刻ずれは0。
生CSV・stderr・終了コード・構築コマンドをgzipとJSONで保持した。
Rustコアとexamplesのソースは変更しておらず、前回の276テスト通過記録を参照する。
この実験のためにRustチェックを再実行したとは主張しない。

~~~bash
python3 -B scripts/run_native_ppc_seed_comparison.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --native-build /workspace/.cloud-setup/solver-reference-20261021 \
  --output-dir /tmp/new-native-seed-comparison

python3 -B docs/comparison/2026-10-21-native-default-seed/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
~~~

次はRustにもSPPによる移動時の初期化を明示的な選択肢として接続し、
同じPPC入力で失敗・再捕捉・誤差への影響を確認する。
同時に、残る大きな外れと移動FIX・単周波制限は別に解決する必要がある。
確認済みチェックリストは18/24（75%）。
