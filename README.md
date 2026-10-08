# gnss-rust

[libgnss++](https://github.com/rsasaki0109/gnssplusplus-library) のRust移植。
移植元は `develop` のコミット
`72f3b7c2c4c088dc575286a6f0edf4e407bc499a` に固定しています。

最終目標は **RTK・PPPが元のlibgnss++と同等の精度・解の状態・処理性能を持つこと**
です。[GOALS.md](GOALS.md)に比較基準と未完了項目を記録しています。
現時点ではコア機能、RINEX 3読み込み、放送暦、初期SPP、LAMBDA候補探索、
初期GPS L1静止FLOATと、全衛星の整数候補を検証する初期FIXを実装しています。
PPPの土台としてSP3/CLK精密プロダクトの読み込み・補間にも対応しています。
送信時刻の精密衛星状態と、明示的な信号補正を使う初期コードSPPも実装しています。
二周波の電離層フリーコード・位相観測モデルも実装しています。
位相アークの継続・再初期化と、採用した観測だけを履歴へ反映する状態管理に対応しています。
これらを接続した初期GPS静止二周波PPP FLOATと、天頂対流圏遅延の推定を実装しています。
Bias-SINEX 1.00の衛星OSBを読み、明示的な時計・位相基準の宣言でPPPへ接続できます。
太陽方向による公称姿勢の位相wind-up補正も、PPP設定で明示的に有効化できます。
絶対ANTEX 1.4の受信機PCO/NOAZI-PCVも、校正と設置オフセットを指定してPPPへ接続できます。
GPSの衛星IF-PCOも公称姿勢で回転し、COM暦からAPCでの測定モデルへ接続できます。
Love数による近似固体地球潮汐を明示的に選び、受信機の基準座標と瞬時座標を区別できます。
IERS2010固体地球潮汐のStep-1＋Step-2も、天体座標・TT/UTを明示入力する部品として実装しています。
GPS→UTC/TAI/TT/UT1とERA・極運動も、期限付きうるう秒表と明示EOPで評価できます。
RTK・PPP双方の同等性能という最終目標は未達です。

## 実装済み

- WGS-84、光速、自転角速度、搬送波周波数・波長、電離層フリー係数
- 衛星IDの表示・解析・順序、信号種別、GLONASS FDMA周波数
- GPS週・週内秒の正規化、符号付き加算、時刻差、GPS尺度での`SystemTime`変換
- ECEF ↔ 測地座標、ECEF差分 ↔ ENU、Sagnac補正付き幾何距離
- 欠測を`Option`で表す測位結果の基本型と有効性・FIX状態の判定
- GPS/Galileo/QZSS/BeiDou放送暦、相対論補正、衛星位置・速度・時計
- RINEX 3観測のストリーム読み込みとKepler航法メッセージ読み込み
- 追尾コードごとのコード・位相・Doppler・SNR、LLI、GLONASS周波数チャネルの保持
- カレンダーからGPS時刻への変換、BDTオフセット、通常のUTC日時のうるう秒補正
- Klobuchar電離層モデル、Saastamoinen対流圏モデル
- GPS/QZSS L1コードSPP、送信時刻・TGD補正、標高角重み付きHouseholder QR
- RINEXファイルからCSV測位結果を出力するCLI
- LAMBDA縮約・MLAMBDA上位整数候補探索、残差/ratioと縮約行列の診断
- 相関した雑音共分散を使うカルマン測定更新、NISゲート、Joseph共分散更新
- 初期GPS L1 C/A静止基線FLOAT、単差アンビギュイティ保持、二重差変換とLAMBDA接続
- 整数候補による位置・共分散の条件付け、残差・履歴検証、連続確認後の初期FIX出力
- SP3-c/d P精密暦・RINEX 3 AS精密時計、最大10点補間、速度・時計ドリフトと相対論補正
- Niellマッピング、初期GPS静止IF PPP FLOAT、実数アンビギュイティ・天頂全遅延状態
- Bias-SINEX 1.00 OSB/DSB読込、期間・時刻系・単位・不確かさ保持、衛星OSBのPPP補正変換
- 公称yaw姿勢のWu位相wind-up、元の近似太陽位置、採用アークごとのcycle継続値
- 元の代替経路に対応するLove数固体地球潮汐・近似月位置、PPPの瞬時受信機座標への接続
- IERS2010 Dehant Step-1＋Step-2の明示天体入力部品、各補正項と合計変位の診断
- GPS時刻尺度変換、うるう秒のUTC準MJD、ERA・極運動・明示歳差章動行列との合成
- 精密衛星の送信時刻評価、信号/時計基準ごとのコード補正、精密コードSPPの合成検証
- 二周波コード/位相の電離層フリー組合せ、信号補正・相関雑音・GF/MW診断
- 位相アークIDと採用回数、GF/MWスリップ候補、欠測・イベント・更新棄却後の再初期化

標準ライブラリのみを使用しています。角度はラジアン、距離はメートル、
周波数はHz、時間は秒です。UTM変換は未実装です。

## 使用例

```rust
use gnss_rust::{GeodeticCoord, GnssTime};
use gnss_rust::coordinates::{geodetic_to_ecef, ecef_to_geodetic};

fn main() -> Result<(), gnss_rust::Error> {
    let position = GeodeticCoord::new(
        35_f64.to_radians(), 139_f64.to_radians(), 45.0,
    )?;
    let ecef = geodetic_to_ecef(position)?;
    let round_trip = ecef_to_geodetic(ecef)?;
    println!("ECEF: {ecef:?}, height: {} m", round_trip.height);

    let time = GnssTime::new(2300, 604_799.5)?.checked_add_seconds(1.0)?;
    assert_eq!((time.week(), time.tow()), (2301, 0.5));
    Ok(())
}
```

ENU変換の入力は「対象位置 − 原点位置」のECEF差分です。
原点のECEF座標をそのまま渡すと意図した相対位置にはなりません。

## 開発・検証

ツールチェーンはRust 1.99.0に固定しています。このクラウド環境では、
各シェルで最初に次を実行してください。

```bash
source /workspace/.cloud-setup/activate.sh
cd /workspace/gnss_rust
cargo test --locked --offline
cargo fmt --check
cargo clippy --locked --offline --all-targets -- -D warnings
cargo build --locked --offline --release
```

通常のRust環境では最初の`source`は不要です。
テストは極・赤道・両半球・衛星高度、週跨ぎ・負の時刻・小数秒、
非有限値・範囲超過、欠測結果を扱います。元のC++ヘッダーから生成した
8ケースの座標基準値、20ケースの放送暦状態、4ケースの大気モデルとも比較します。
[基準値の再生成手順](tests/fixtures/README.md)
を参照してください。Rustのビルド・テストにC++やEigenは不要です。

## RINEXからSPPを実行する

大気モデル込みの合成観測を使う例です。

```bash
cargo run --release --locked --offline --bin gnss-rust -- spp \
  --obs tests/fixtures/synthetic_spp_atmosphere.obs \
  --nav tests/fixtures/synthetic_spp.nav
```

標準出力はGPS週・週内秒、ECEF、時計バイアス、衛星数、状態、残差RMSなどのCSVです。
`--output new-results.csv`で新規ファイルに保存できます。既存ファイルは上書きしません。
大気遅延を含まない`synthetic_spp.obs`では`--no-atmosphere`を指定します。
合成入力は9衛星・8エポックで、既知の真値から各ECEF成分3 mm以内、
時計バイアス差1e-11秒以内、残差RMS1 mm未満をテストします。
これはノイズのない合成データの回帰検証で、実測データでの精度保証ではありません。

SPPは各衛星のC1Cを優先し、なければ追尾コードの辞書順でL1コードを選びます。
GPSとQZSSは共通時計を推定します。未対応信号、欠測コード、未収録/古い/不健康な暦、
標高角マスクによる除外は結果に記録します。階数不足や衛星数不足はエラーです。
CLIは測位エラー時にエポックを示して非ゼロ終了し、成功した実行だけ完了件数を表示します。
途中失敗時は標準出力や新規出力ファイルに途中までのCSVが残る場合があります。
flag 6のサイクルスリップ記録は測位せず、スキップ件数を報告します。

## 現在の対応範囲

- 観測はRINEX 3。通常・電源異常・サイクルスリップのフラグを保持し、
  ヘッダー変更イベントと観測/ヘッダーの継続行を扱います。
- 航法はGPS・Galileo・QZSS・BeiDouのKeplerメッセージ。BeiDou GEOの座標枠も扱います。
  GLONASS/SBAS航法、RINEX 2/4、圧縮RINEXは明示的に未対応です。
- 非単位スケール係数と非ゼロ位相シフトは黙って無視せず拒否します。
  受信機時計オフセットはメタデータとして保持し、時刻・コードへ追加適用しません。
- SBAS/QZSSの衛星番号再マッピング、Galileo SSRに応じたI/NAV/F/NAV選択、
  信号別バイアス選択は未完了です。RINEXの衛星番号はそのまま保持します。
- 通常のUTC日時は2017-01-01までの歴史的うるう秒表を使います。
  以降は最後の既知オフセット18秒を使い、将来の変更を予測しません。
  うるう秒そのものの`second == 60`は未対応です。
- SPPはGPS/QZSS L1。高度な外れ値除去、RAIM/FDE、二周波、速度、
  マルチGNSS時計、SSR、元のSPPProcessor全体の処理は未移植です。
  精密コードSPPは下記の別API・例で実行します。
- RTKのFIX保持/移動体/部分AR/FFRT、PPPのマルチGNSS/移動体/AR、衛星PCV/IERS2010の自動入力・接続、海洋荷重・極潮汐等は今後の実装対象です。

## 初期RTK FLOATを実行する

```bash
cargo run --release --locked --offline --bin gnss-rust -- rtk-float \
  --rover tests/fixtures/synthetic_rtk_rover.obs \
  --base tests/fixtures/synthetic_rtk_base.obs \
  --nav tests/fixtures/synthetic_rtk.nav \
  --base-position '-3947481.0649388283,3431492.9375319984,3637892.7203177349'
```

上の座標は合成基準局の真値です。`tests/fixtures/synthetic_rtk_truth.csv`の
最初のデータ行の3列にも記録しています。CLIには`--base-position X,Y,Z`の2引数で
指定します。実データの基準局には測量されたWGS-84 ECEF
座標が必要です。移動局の初期座標はRINEXヘッダーから読み、`--rover-position X,Y,Z`で
指定できます。出力は位置・衛星数・FLOAT・基準衛星・残差・NIS・再初期化数・候補ratioのCSVです。
`--output NEWCSV`はSPPと同様に新規ファイルだけを作成します。

現在の処理は**GPSのC1C/L1C、静止基線、既定の基線長上限10 km**に限定しています。
二重差は`(rover − base)reference − (rover − base)satellite`です。単差アンビギュイティを
保持するため、基準衛星変更の際に状態全体を初期化しません。雑音共分散には共有基準衛星の
相関を含めます。LAMBDAに渡す二重差状態と位置との交差共分散はAPIで取得できます。

時刻差は1e-7秒以内、GPS時刻は単調増加、両ファイルのエポック列は一致する必要があります。
補間や無言のエポック間引きはしません。非ゼロ受信機時計オフセットを含む入力は未対応です。
LLIのloss-of-lock、欠測後の再捕捉、power-failure、既定120秒超の間隔でアンビギュイティを
再初期化し、half-cycleは除外します。flag 6の同期イベント対はスキップして継続性を切ります。
更新が失敗すると直前の数値状態を保持し、次の有効対でアンビギュイティを再初期化します。
CLIは不一致・測位失敗時に非ゼロ終了し、途中までのCSVが残る場合があります。

`rtk-float`のratioは診断値で、**全エポックがFLOATのまま**です。候補探索が失敗した場合はratioを空欄にし、
理由と診断欠測数を標準エラーへ報告します。FIX保持、部分AR/FFRT、LLIなしの
スリップ検出、移動体、多周波・マルチGNSS、電離層推定、外れ値の個別除外は未移植です。
対流圏はSaastamoinen、電離層は推定せず短基線で相殺すると仮定しています。
基線長の上限だけで大気誤差の影響が消えることは保証しません。時計は二重差で消去され、
相対測位では受信機時計バイアスを推定しません。

合成データ32エポックでは3次元誤差3 mm未満をテストし、コード雑音振幅5 cm・位相雑音
振幅0.5 mmを含む32エポックでは15 cm未満をテストします。実データ精度、FIX率や
元のRTKProcessor全体との同等性を示す結果ではありません。

## 初期RTK FIXを実行する

上のコマンドの`rtk-float`を`rtk`に変更すると、整数候補の検証を有効にします。
入力条件は同じです。全組の二重差整数で位置と共分散を条件付けし、同じ送信時刻・
観測・衛星選択でコードと位相の残差を再計算します。以下の既定条件をすべて満たし、
同じ整数候補を連続して確認したときだけFIXを出力します。

| 判定 | 既定条件 |
| --- | --- |
| 観測・候補 | GPS L1、二重差4組以上、全衛星5エポック以上のロック |
| LAMBDA | 最良2候補のratio ≥3、探索完了 |
| 条件付き位置 | FLOATからの変位 ≤1 m、各ECEF成分の標準偏差 ≤5 cm、基線 ≤10 km |
| 位相残差 | RMS ≤1 cm、最大絶対値 ≤3 cm |
| コード残差 | RMS ≤1 m |
| 位置に依存しないコード・位相整合 | 全組で `DD_PR − (DD_CP − λN)` の絶対値 ≤0.5 m |
| 残差の共分散重み付き二乗和 | 観測行あたり ≤9 |
| 静止位置の履歴 | 120秒以内の直前FIXとの変位 ≤10 cm |
| 連続確認 | 同じ衛星組の基準衛星に依存しない整数ラベルを2エポック以上で確認 |

これらは初期Rust実装の保守的な設定です。元RTKProcessorの全設定・FIX状態遷移を
再現するものではなく、誤FIX率の保証でもありません。共分散の対角値を下限で補正せず、
結合共分散が正定値でなければ棄却します。最良LAMBDA残差が厳密にゼロの場合は、
元のratioゼロ規約を維持するためFIXを宣言しません。APIの`RtkFixConfig`で設定できます。

各エポックの検証に失敗した場合はFLOATを出力し、連続確認をやり直します。
LLI・再捕捉・電源異常・長い間隔ではロックを再計数します。FIXの位置や整数をFLOAT
フィルターへ戻さず、整数の保持も行いません。`RtkFixPolicy::evaluate`はFLOATの
不変な観測スナップショットを検証し、`validate_candidate`は任意整数の診断に使えます。

CSVの位置・コード/位相RMSは出力状態に対応します。`float_update_nis`はFLOAT更新の
NISです。追加の`fix_decision`は棄却・確認・採用理由、`candidate_pairs`は探索した組数、
`confirmations`は連続確認数、`candidate_postfit_nis_per_row`は候補の残差検証値です。
ロック不足など探索前のエポックはratioと候補残差が空欄です。

同じ32エポックの合成入力では、初めの5エポックがFLOAT、残る27エポックがFIXに
なります。最大3次元FIX誤差は雑音なし約0.032 mm、雑音あり約0.86 mmでした。
全FIXの整数が既知真値に一致すること、誤整数や高ratioでも位相残差の大きい候補の棄却、
スリップ後のFLOAT・再確認、基準衛星変更、失敗後の復帰をテストします。
実データ・元RTKProcessor全体・処理速度/メモリの比較は未実施です。

## LAMBDA整数候補探索

```rust
use gnss_rust::lambda::{search, LambdaConfig};

fn main() -> Result<(), gnss_rust::lambda::LambdaError> {
    let result = search(
        &[1.4, 2.4],
        &[vec![1.0, 0.99], vec![0.99, 1.0]],
        LambdaConfig::default(),
    )?;
    println!("best: {:?}, ratio: {:?}", result.candidates[0], result.ratio());
    Ok(())
}
```

入力はcycle単位のfloat ambiguityとcycle²単位の共分散です。共分散は行優先の
対称正定値行列です。候補は整数ベクトルとMahalanobis残差の小さい順に返し、
縮約後の条件付き分散と変換行列も保持します。独立な四捨五入では誤る相関を扱います。
これだけで位置やFIXを出すことはありません。

元実装と同じ`floor(x + 0.5)`、同点候補の順序、最良残差ゼロ時のratioゼロを保持します。
逆変換はチェック付き整数演算に変更しています。1〜128変数・1〜64候補に対応し、
縮約・探索の既定上限はそれぞれ10,000回です。整数の絶対値が2⁵²以上になる演算、
非有限値、不正な共分散、探索上限はエラーです。途中の候補を成功扱いしません。
元実装の縮約比較定数1e-6と初期探索距離1e99も保持しています。
FFRT、bootstrapped success rate、部分AR、FIX保持は今後の実装対象です。
11ケース・34候補の元C++比較と、相関共分散100ケースの総当たり照合をテストします。

## SP3精密暦・CLK精密時計

`precise::read_sp3`と`read_clk`で読み、`PreciseProducts::interpolate`でGPS時刻の
位置・速度・公開時計を求めます。合成ファイルを使った実行例です。

```bash
cargo run --release --locked --offline --example precise -- \
  tests/fixtures/synthetic_precise.sp3 tests/fixtures/synthetic_precise.clk
```

例は最初の暦間隔の中点（GPS週2300、週内秒345750）で全衛星を評価し、
ECEF・速度・時計の出所・相対論補正をCSVに出力します。CLK引数を省くとSP3時計を使います。
この例は衛星状態の評価です。コード測位への接続は次節、位相PPPへの接続は
下記の初期PPP API・例で行います。RTKへの精密プロダクト接続は未実装です。

- SP3-c/dのPレコード、RINEX CLK 3.00〜3.04のASレコードに対応します。
  CLKは1〜6値、継続行、FortranのD指数を扱い、時計sigma・drift・rateも保持します。
  SP3のV/EP/EV、旧SP3、ヘッダーなしSP3、非AS時計、圧縮ファイルは未対応です。
- SP3のkm・microsecondをm・sへ変換します。CLKはsです。ゼロ時計は有効で、
  SP3時計の欠測値は`None`です。ゼロ座標成分は保守的に欠測として扱います。
  宣言衛星数・エポック数・間隔・重複・不正値を検査し、不完全なSP3を成功扱いしません。
- GPS/GAL/QZS/BDT/UTCの時刻系をGPS時刻へ変換します。CLKの時刻系省略は
  仕様のGPS既定値として記録します。TAI/GLO・うるう秒そのものは未対応です。
  UTCの将来うるう秒は既存の時刻モデルと同じく予測しません。
- 座標系、COM/APC（コメントによる位相中心宣言）、精度指数・基数、品質フラグ、
  生ヘッダーを保持します。座標枠変換や精度指数からの重み生成は行いません。
- 元C++と同じ最大10点のNeville多項式補間、±1秒の中央差分で速度・時計ドリフトを
  求めます。最低2点が必要です。既定では外挿を行わず、900秒を超えるサンプル間隔、
  欠測、SP3の軌道manoeuvre・時計イベントをまたぐ補間も行いません。
  `InterpolationConfig`で外挿を最大900秒まで明示できますが、内部欠測は橋渡ししません。
  予測データは既定で除外し、許可時には結果の補間点と予測利用を記録します。
  SP3品質フラグはSP3系列に適用します。CLKに旗のない時計ジャンプの自動検出は未実装です。
- 衛星にCLKがあれば、その衛星の時計系列全体にCLKを選びます。CLKの範囲外や欠測で
  SP3時計へ切り替えず、`clock: None`と理由・選択した出所を返します。CLKにない衛星はSP3時計を使います。
  元C++のSP3/CLK混合系列との差をテストで記録しています。
- 時計は公開値をそのまま返します。`precise_clock_relativistic_correction`は
  `−2(r·v)/c²`を別に計算します。将来の測位では送信時刻で一度だけ時計に加算し、
  アンテナ・信号バイアス等の必要な補正も別に扱う必要があります。

C++比較64ケースでは位置・速度・相対論補正を比較し、同じ時計系列の55ケースは
時計・ドリフトも比較します。残る9ケースは元C++の混合時計との差を記録するケースです。
2点の解析解、週跨ぎ、小数秒、欠測、品質境界、時刻系、CLK継続行も検証します。
合成多項式データの部品比較で、実プロダクトの精度やPPPソルバーの同等性ではありません。

## 送信時刻評価と精密コードSPP

`precise_transmit::evaluate_precise_transmit`は受信GPS時刻と受信機位置の推定値から
送信時刻の精密衛星位置・時計を求めます。既定の`ReceptionTaylor`は元C++の精密経路と
同じく、受信時刻の状態に対して`p − vτ`、`clock − driftτ`を計算します。
`IteratedEmission`を指定すると送信時刻で再補間し、幾何学的な光行時間を反復します。
両方式の差をテストに残し、一次近似の結果を厳密な再補間と同一とは扱いません。

τは元C++と同じECEFのユークリッド距離/cです。Sagnacはその後の測定モデルで一度だけ
加算します。公開時計と周期相対論項は結果に別々に保持し、`corrected_clock_bias_s`は
両者の和です。放送暦時計、TGD、アンテナ、コード/位相バイアスは自動適用しません。
受信・送信の両端点で精密暦と選択時計が使える必要があり、品質境界をまたぐ投影も拒否します。
元C++のファイル端での投影とは採用条件が異なります。光行時間上限は既定1秒、
再補間反復は既定8回・収束閾値1e-11秒で、未収束はエラーです。

`precise_spp::solve_epoch_precise`でGPS/QZSS L1のコード位置と受信機時計を推定できます。
以下の合成例では、9衛星・8エポックのRINEX観測、精密暦/時計、既知のコード補正を使います。

```bash
cargo run --release --locked --offline --example precise_spp -- \
  tests/fixtures/synthetic_precise_spp_taylor.obs \
  tests/fixtures/synthetic_precise_spp.sp3 \
  tests/fixtures/synthetic_precise_spp.clk \
  tests/fixtures/synthetic_precise_spp_bias.csv
```

再補間方式を試す場合は観測を`synthetic_precise_spp_emission.obs`へ変更し、末尾に
`emission`を追加します。この例はIGS20/APC・電離層なし・Saastamoinen対流圏ありの
合成モデル専用です。通常の`spp` CLIは従来の放送暦経路を使います。

各信号のコードには`P_corrected = P_raw + correction_m`という符号で補正を指定します。
`PreciseCodeBiases`には衛星、追尾コード、時計の出所、時計基準名、受信GPS時刻での
有効期間（両端を含む）を宣言します。ゼロ補正も明示的なエントリが必要です。
欠測・期限切れ・時計の出所/基準名の不一致は棄却し、採用した補正と送信時刻状態を結果に
保持します。同じ衛星・コード・時計の出所で重なる期間は拒否します。
時計基準名は呼出側の宣言で、ファイルからの自動検証ではありません。

例の補正CSVはテスト用の簡単な台帳です。下記のBias-SINEX読込は初期PPPへ接続します。
時計基準の自動変換は未実装です。実プロダクトのC1Cへ放送暦TGDをそのまま流用してはいけません。
`PreciseSppConfig`は期待する座標系・COM/APCを明示してヘッダーと照合し、
座標枠やアンテナの変換は行いません。電離層モデルを有効にする場合は係数を明示します。
非ゼロの未適用受信機時計オフセットは拒否します。

C++補間・幾何・対流圏部品で生成した合成観測では、両方式とも8エポックすべてがSPPで
既知位置・時計を回復します。最大3次元誤差は一次近似約0.0107 mm、再補間約0.0098 mm
でした。送信時刻80ケースの部品比較、相対論の符号と一度だけの適用、時計欠測、品質境界、
補正の不一致、週跨ぎ、初期座標なしも検証します。元SPPProcessor全体の互換性や実測精度を
示す結果ではありません。二周波・位相・ambiguity状態を持つ初期PPPは下記の別APIです。

## 二周波コード・位相観測モデル

`dual_frequency::form_ionosphere_free`はRINEX観測エポックから、指定した衛星・
2つの追尾コードの電離層フリー（IF）コードと位相を生成します。信号の優先順位や
単周波への切替は行いません。GPS/QZSS L1・L2・L5、Galileo E1・E5a・E5b・E6、
BeiDou B1I・B2I・B3I・B1C・B2a、GLONASS FDMA L1・L2の明示的な追尾コードの
部分集合に対応します。対応する文字は`tracking_frequency_hz`に記載しています。
GLONASSは両観測に同じ有効チャネルが必要です。この周波数対応はGLONASS測位や
マルチGNSS受信機時計の実装を意味しません。

係数は`a = f1²/(f1²−f2²)`、`b = −f2²/(f1²−f2²)`です。
コードは`a P1 + b P2`、位相はcycleから各波長でメートルへ換算して
`a λ1 L1 + b λ2 L2`とします。一次の分散性電離層項は相殺しますが、位相には
`a λ1 N1 + b λ2 N2`という実数のメートル単位アンビギュイティが残ります。
これを整数cycleとしてLAMBDAへ渡してはいけません。高次電離層は消去しません。

`SignalCorrections`で両信号のコード加算補正を宣言し、位相を使う場合は位相加算補正と
位相基準名も宣言します。ゼロも明示します。コード・位相とも補正単位はメートルで、
IF組合せの前に加算します。衛星・追尾コード・時計出所・時計基準・位相基準と
受信GPS時刻の有効期間（両端を含む）を照合します。欠測コード、コード補正の欠測、
期限切れ、基準不一致、不正な観測や重複はエラーです。同じ信号・時計出所の期間重複も
拒否します。これらは呼出側の宣言です。Bias-SINEX OSBからの補正変換は下記のAPIで行い、
実際の精密時計の出所はPPP側で確認します。時計/位相基準の互換性の自動判定・変換、
TGD・アンテナ・wind-upの自動適用は行いません。
PPP側の明示的なwind-up処理は下記の別設定です。

片方の位相欠測、半サイクル、位相補正未宣言ではIFコードを保持し、位相を`None`と
理由で返します。両方の位相が使える場合、LLIのロック喪失またはepoch flag 1で
`reset_required`を返します。生観測・LLI・適用補正も保持します。GFは未補正の
`λ1 L1 − λ2 L2`、MWは未補正のcycle位相とコードから計算し、補正変更による
人工的なジャンプを混ぜません。元PPPの通常経路と同様、GLONASSではMW診断を返しません。
このAPIは状態を持ちません。継続管理とGF/MW判定は次節の`PhaseArcTracker`で行います。
flag 6と非ゼロの未適用受信機時計オフセットは、この測定ビルダーでは入力として拒否します。

雑音は各周波数のコード/位相それぞれにメートル²の2×2正定値共分散を指定します。
周波数間の相関を含む`wᵀ Q w`を伝播し、コードと位相の間の交差共分散は表現しません。
分散を一律の下限値で置き換えず、不正な共分散を拒否します。

```bash
cargo run --release --locked --offline --example dual_frequency
```

例は同梱の合成RINEXからG01のC1C/L1C・C2W/L2Wを読み、ゼロ補正を明示して4エポックの
IFコード・位相・位相−コード差・GF/MWを出力します。測位結果ではありません。
元C++の未変更PPP係数/MW関数と観測式による5衛星系・12信号対・48ケースを比較し、
同じ観測をRINEXに丸めて読み戻した48件も検証します。電離層消去、残るメートル単位
アンビギュイティ、補正の符号、相関雑音、欠測・LLI・期間・週跨ぎもテストします。
このモデルは下記の初期PPPフィルターへ接続しています。精密コードSPP/RTKは別経路です。

## 位相アークの継続管理

`phase_arc::PhaseArcTracker`は、各衛星の選択した二周波ペアについて、アンビギュイティを
継続できる区間をアークIDと採用回数で管理します。位置・時計・対流圏・アンビギュイティ値
そのものは推定しません。準備と採用を次の順に分けます。

1. `prepare_epoch`へ観測、信号補正、衛星ごとの`DualFrequencyConfig`を渡します。
   IF測定、提案アークID、再初期化理由、比較可能なGF/MW差を返します。
2. 呼出側の推定器がコード・位相の更新を検証します。
3. `finish_epoch`へ実際に位相更新を採用した衛星の集合を渡します。コードだけの採用は
   含めません。全更新失敗時は空集合を渡します。

準備しただけではアーク・ロック数・GF/MW履歴を更新しません。提案は内部に保持し、
確定前の次エポック処理を拒否します。不適格な衛星の採用指定は訂正できます。
未採用の衛星は継続を切り、次の有効位相は新しいアークから開始します。既存数値フィルター
状態の維持やアーク変更時の状態/共分散の再初期化は呼出側の仕事です。
終了済みエポックは失敗時も時刻を進め、同じエポックの重複採用を拒否します。

GFは元C++通常PPP経路の`max(configured_threshold, 0.5)` m、MWは
`max(configured_threshold * 100, 10)` mを用い、絶対差が閾値を**超えた**場合に
新しいアークを提案します。既定の設定値は0.05 mなので、実効閾値はGF 0.5 m・MW 10 mです。
等号は継続します。両端にMWがなければMWで判定せず、GLONASSはGFだけを使います。
比較は未補正の診断量で行い、既知の欠測・信号変更・間隔・LLI等をまたぐ差を使いません。
MWはコード異常でも変化するため、判定理由はジャンプ候補であってスリップの真値ではありません。

片方の位相欠測、half-cycle、衛星欠測、選択対象からの除外、測定生成の棄却、
更新棄却は継続を切ります。LLIのロック喪失、電源異常、既定120秒を**超える**間隔、
追尾コード/周波数/GLONASSチャネル・時計出所/基準・位相基準の変更も再初期化します。
この初期政策はコード/位相の加算補正値の変更にも保守的に再初期化し、補正の有効期間延長や
雑音設定だけの変更では継続します。wind-upは下記のPPPで校正値とは別に適用し、
滑らかな変化だけではアークを切りません。汎用的な動的アンテナ補正の継続管理は今後の
対象です。flag 6はイベント専用として、記録にない衛星も含め全管理対象の継続を切ります。
不正な前進エポックは既存アークを切ってエラーにします。重複・逆順エポックは状態を変えず
拒否します。観測列の外で継続が失われた場合は`invalidate_all`を使います。

```bash
cargo run --release --locked --offline --example phase_arc
```

例は2衛星・16エポックの合成RINEXを読み、提案/採用アークID、採用回数、理由をCSVに
出力します。8番目のエポックではG01の位相更新を意図的に棄却し、G02だけを継続します。
合計26件の位相採用と13個の採用アークを検証しています。これは継続管理の実行例で、
PPP推定や測位結果ではありません。C++のGF閾値ヘルパーと通常PPPの比較式を使う48件の
閾値参照は、正負・等号・直上・MW欠測を含みます。元PPPProcessor全体のスリップ政策、
SSR discontinuity counter、CLAS/MADOCA専用政策、実データの誤検出率は未検証です。

## 初期GPS静止PPP FLOAT

```bash
cargo run --release --locked --offline --example ppp
cargo run --release --locked --offline --example ppp -- noisy
```

同梱の合成SP3/CLK・GPS C1C/L1CとC2W/L2W・Bias-SINEX OSBを使い、各64エポックの
ECEF・受信機時計・天頂全遅延・コード/位相残差・使用数・NIS・アーク再初期化数を
CSVへ出力します。この例の引数は合成データの選択だけで、汎用実データCLIではありません。

`ppp::PppFilter::new(PppConfig)`には独立した初期座標・時計・天頂遅延、SP3座標系名と
COM/APC宣言を必須で渡します。座標系の自動変換は行いません。COM入力には下記の
`satellite_antenna`を必須にし、APC入力には`None`を指定します。
受信機には下記の`receiver_reference`と`solid_earth_tide`を指定します。衛星COM/APCの
宣言とは独立しており、潮汐を使わない従来経路は`InstantaneousMarker`/`None`です。
`process_epoch`へ観測・精密プロダクト・`SignalCorrections`・衛星ごとの
`DualFrequencyConfig`を渡します。GPS PRN 1〜32の単一受信機時計に限定し、他衛星系の
選択を明示的に拒否します。信号/時計/位相基準・有効期間の宣言に加え、実際の精密時計の
出所との一致を確認します。選択したCLKが欠測でもSP3時計や放送暦へ切り替えません。

状態は`[ECEF XYZ, 時計(m), 天頂全遅延(m), 衛星別IFアンビギュイティ(m)]`です。
測定式は`rho + 時計 - c*衛星時計 + Niell*天頂全遅延`、位相には実数アンビギュイティを
加えます。周期相対論とSagnacはそれぞれ一度だけ適用し、放送暦TGDは加えません。
一次電離層はIFで消去します。天頂全遅延をNiellのhydrostatic mappingで写す初期経路で、
乾燥/湿潤の分離・勾配推定は未実装です。日番号は元の現在の時刻変換と同じGPS尺度の
日付を使います。XYZの設計行には`-LOS`に一次Sagnac微分を加え、マッピングと精密衛星
位置と潮汐変位の受信機座標微分は省略します。元通常PPPの`-LOS`だけの行とはこの点が異なります。

既定の初期分散はXYZ各`1e6`、時計`1e8`、天頂全遅延`100`、新規アンビギュイティ
`1e14` m²です。アンビギュイティは独立なゼロ平均・広い事前分布から開始し、位相
カウンターが約−2,000万mずれる初期ゲージも検証しています。同じエポックの位相やSPP
結果を独立な事前情報として再利用しません。静止XYZの既定プロセス雑音はゼロ、時計は
`100`、天頂全遅延とアンビギュイティは各`1e-8` m²/sで、最後の数値更新からの経過時間を
一度だけ積算します。各反復は同じ予測事前分布から再線形化し、観測情報を重複加算しません。
周波数間の相関雑音はIFへ伝播しますが、コードと位相の間・衛星間の測定雑音は独立とします。

既定は標高角10°、独立な幾何を持つ5コード以上・4位相以上、最大8反復、収束変化
`1e-4` m、NIS/測定行数≤25、更新後絶対残差はコード≤2 m・位相≤0.05 mです。
反復中の採用衛星は固定します。コードだけの解をPPPとして出力せず、成功状態は
`PppFloat`です。整数AR/FIX・速度・外れ値衛星の逐次除去は未実装です。

アークIDが継続する衛星は状態と相関共分散を保持し、欠測・再初期化衛星の古い状態は
残る状態の共分散部分行列を保って除去します。新規状態は独立に追加します。数値更新が
棄却された場合は平均・共分散・数値更新時刻をそのまま保持し、位相継続を切ります。
`state()`はその最後の成功時点のスナップショットなので、時刻とアーク状態も確認してください。
イベントは継続管理だけを更新します。重複・逆順の入力は状態を変えず拒否します。
成功結果には補正・精密送信状態・残差・実数アンビギュイティ・採用アークを残します。

元C++のNiell 96件、精密幾何/時計/測定式54件、全64エポック×9衛星の標高角採用576件を
比較しています。合成データの全64エポックがPPP FLOATとなり、最後の3次元位置誤差は
雑音なし約0.009 mm、コード振幅5 cm・位相振幅1 mmの決定的雑音あり約2.66 mmでした。
最初の56エポックは9衛星、最後の8エポックは標高角マスクで8衛星を使用します。
棄却後の数値状態維持、通知あり/なしスリップ、位相/時計欠測、基準/補正変更、イベント、
悪い幾何・非収束・重複入力も検証します。これは部品比較と合成真値への回帰結果です。
元PPPProcessor全体の状態・共分散・収束政策、実データ精度、速度/メモリの同等性は未確認です。
衛星PCV・潮汐、実際のyaw姿勢/食への対応、差分バイアスの基準変換などは今後の対象です。

## Bias-SINEX OSBをPPPへ接続する

`bias::read_bias_sinex(reader, BiasReadOptions)`はBias-SINEX 1.00の`BIAS/SOLUTION`を
読みます。OSBとDSB、SVN/PRN、追尾コード、受信機局名、有効期間、単位、値とsigmaを
保持します。空のOSB OBS2列とSVNなしの短縮行、OBS2に同じ観測名を繰り返す元実装の
形式を扱います。記述ブロックの時刻系・時計参照観測なども保持し、衛星系ごとの同名の
記述行を上書きしません。OSBの異なる2観測名は意味が曖昧なので補正変換時に拒否します。

`TIME_SYSTEM`はG/U/E/C/Jに対応します。記載がない場合は呼出側の時刻系指定が必須で、
指定との矛盾も拒否します。UTC/BDTはGPSへ変換し、年通日と秒の範囲を確認します。
`0000:000:00000`は明示的な開いた端点です。期間は**開始を含み終了を含まない**ので、
隣接する期間の境界では新しい行だけを採用します。重複・期間重複・不正な日付・非有限値・
負sigma・未対応の単位やbias slope・途中で切れたブロック/終了印をエラーにします。
IONEX補助DCB、受信機だけの行、他バージョン・圧縮ファイルは未対応です。

単位はnsとm（meter/metersも可）、位相OSBのみcycを扱います。nsは`c*1e-9`、cycは
正確な追尾コードの波長でmへ変換し、GLONASS FDMAのcycにはチャネル指定が必要です。
元PPPの通常コードOSBと同じ減算符号で、コード/位相の加算補正へ変換します。
位相OSBの適用とcyc換算は今回のRust経路で、元PPPProcessor全体との位相比較ではありません。
DSBだけでは絶対バイアスの基準を決められないため保持だけにとどめ、OSBの代用にはしません。
相対バイアス宣言のある製品や受信機局付き行も衛星絶対補正として採用しません。

`BiasProducts::corrections_at(time, &OsbBinding, &requests)`に製品識別名、互換性を確認した
時計出所/基準名・位相基準名、衛星と追尾コードを渡します。製品の文字列から時計/位相
基準の互換性を自動認定しません。返る`BiasCorrectionSnapshot`は**その受信時刻だけ**の
`SignalCorrections`と、元行・sigma・製品名の来歴を持ちます。PPPの`process_epoch`へ
`snapshot.corrections`を渡し、各エポックで再生成してください。期間が毎回変わっても、
基準と補正値が同じなら位相アークは継続します。補正値の切替は現在の保守的な政策で
該当衛星のアークを再初期化します。sigmaは来歴に残し、毎エポック独立な雑音として
自動加算しません。バイアス誤差の時間相関を持つ推定モデルは未実装です。

要求したコードOSBの欠測/期限切れはエラーです。位相OSBの欠測/期限切れ・位相基準未宣言は
理由を来歴に残し、コードだけを保持します。ゼロ位相補正で埋めません。補正生成で失敗して
PPP処理前にエポックを飛ばす場合、呼出側で`PppFilter::invalidate_phase()`を呼び、継続を
切ってください。使用例もこの経路を使います。衛星除外の判断は呼出側で行います。

元C++の未変更DCBProducts読込39行と、通常PPPの単位変換式・符号を比較しています。
実装元は空OBS2のOSBを取りこぼすため、C++比較用にはOBS2を繰り返した別ファイルを使い、
空欄の標準形式は別の合成ファイルで検証します。同じ36 OSBを使ったPPPの全128エポックで
従来の補正台帳経路と数値状態・共分散・採用アークを照合し、release使用例のCSVも一致しました。
期間境界、位相補正欠測/復帰、時計基準不一致、外部での補正失敗後の復帰も検証しています。
実配布製品の時計/位相基準や実データ精度・全体性能の同等性はまだ確認していません。

## 公称姿勢の位相wind-up

```bash
cargo run --release --locked --offline --example ppp -- windup
cargo run --release --locked --offline --example ppp -- windup noisy
```

元C++のwind-upを両搬送波へ加えた合成観測を選び、対応する補正を有効にします。
CSVには既存15列に`windup_satellites`と`windup_max_abs_cycles`、下記の受信機補正診断2列を追加しています。
`windup`を省く場合は従来の補正なし合成観測を使い、既存列の結果はそのままです。
実入力APIでは`PppConfig.phase_windup`へ
`Some(windup::PppWindupModel::NominalYawApproximateSun)`を渡します。
補正済み観測や補正なしの合成データでは`None`を渡し、二重適用を避けてください。

`windup::nominal_yaw_windup`はWuらの有効ダイポール式を使います。衛星z軸は地心方向、
y軸はzと衛星から太陽への方向の外積、x軸はyとzの外積です。受信機は固定のNorth–West–Up
軸とします。衛星から受信機への視線、符号付き主値±0.5 cycle、前回値に最も近い整数を
使い連続値を返します。整数への丸めは元のC++と同じ、等距離時はゼロから遠い方です。
明示的な姿勢軸を持つ部品API`phase_windup_from_frames`もありますが、PPP経路の受信機
姿勢は固定です。車両姿勢や姿勢観測・衛星ブロック別yaw、食/正午/深夜maneuverは未対応です。

太陽位置は元の`ppp_corrections.cpp`の低次近似とGMST回転を移植しています。
小数秒を含む受信GPS時刻ラベルを直接使う元経路に合わせ、GPS→UTC・UT1/EOP・章動・
極運動の補正は行いません。これは高精度な天体暦ではありません。欠けた太陽位置や退化した
姿勢・ダイポールをエラーにし、元の代替姿勢や前回値の無言の再利用は行いません。

cycle補正はコードへ加えず、OSB適用後の位相へ`-cycles*lambda`を各周波数で適用します。
IFには恒等式`a*lambda1+b*lambda2 = c/(f1+f2)`による減算を一度だけ適用します。
`PppMeasurement.observation`はwind-upと受信機アンテナ補正の適用前、`windup`には太陽位置、前回/主値/継続cycle、反復の整数基準、
周波数別加算補正とIF加算補正を残します。測定の実効位相は
`observation.phase.phase_if_m + windup.phase_if_add_m`に、受信機補正を有効にした場合は
`receiver_antenna.if_add_m`も加えた値です。GF/MWは生観測のままです。
wind-upの受信機座標微分は設計行へ含めません。

反復ごとに現在の位置・精密送信状態から再評価します。継続アークは同じ採用済み前回値を
整数基準に使います。新規アークは予測位置でゼロ基準から一度評価し、その値をエポック内の
反復基準に固定します。半cycleの主値境界で反復ごとに整数が変わることを防ぎます。
反復途中の候補をcacheへ戻しません。全数値検証と位相採用が成功した時だけ、
アークID・時刻・cycleを保存します。スリップ、位相欠測、マスク/除外、イベント、電源異常、
更新棄却・外部無効化で継続を切り、新規アークの整数基準はゼロから開始します。
整数基準の変更は新しい実数アンビギュイティで吸収します。重複/逆順入力は採用済み状態を
変えません。使用位相のwind-upが計算できない場合はエポック更新を棄却し、数値状態を維持します。

実際の元C++太陽ヘルパー12件・wind-up関数360件を比較し、周期3周の正負の回転・
主値分岐跨ぎ・丸め境界・不正な姿勢も検証しました。補正入りの両64エポックでアークは
滑らかに継続し、最後の3次元位置誤差は合成雑音なし約0.012 mm・あり約2.66 mmです。
棄却後のcache未採用、再捕捉、位相欠測、電源異常、疎なイベント、退化姿勢時の更新棄却も
検証しています。元PPPProcessor全体の姿勢/cache政策・実データ精度・処理性能は未比較です。

## 受信機ANTEX PCO/PCV

```bash
cargo run --release --locked --offline --example ppp -- antenna
cargo run --release --locked --offline --example ppp -- antenna noisy
```

`antenna`は受信機PCO/PCVとwind-upを注入した合成観測を選び、両補正を有効にします。
衛星暦は従来のAPCのままです。CSVの末尾に`receiver_antenna_satellites`と
`receiver_antenna_max_abs_if_m`を出力します。引数なし/`windup`だけの既存17列は変わりません。

`antenna::read_antex`は絶対ANTEX 1.4の混在ファイルを読み、完全なアンテナ型名・radome、
serial、SVN、周波数、GPS暦の有効期間、NEUのPCO、NOAZIの全ノードと元行を保持します。
mmをmへ変換し、PCVは元の線形補間と端点クランプを使います。`Antex::receiver`へ
型名・serial・時刻を渡し、期間内の唯一の受信機校正を選びます。空serialの汎用校正も
明示的に指定し、型名/radome/serialの不一致や複数候補から自動で代替しません。
有効期間の両端は含み、省略端は無制限です。相対PCV、方位依存DAZI、周波数RMS、
欠けたPCO/NOAZI、誤ったノード数・周波数数、未完了ブロックはエラーにします。

`ReceiverAntennaModel::new(product_id, calibration, marker_delta_enu_m)`で選択を固定し、
`PppConfig.receiver_antenna = Some(model)`を設定します。設置オフセットはmarker→ARPの
**East/North/Up（m）**、ANTEXのPCOは**North/East/Up（m）**です。GPS 1/2/5の追尾コードを
ANTEX G01/G02/G05へ照合します。衛星側の校正は読込後も受信機として選択しません。
RINEXヘッダーの自動照合や既適用補正の検出は行わず、補正済み観測では`None`を指定します。

markerから衛星へのENU単位視線を`u`とし、各周波数で
`add = u·(marker_delta_enu + PCO_enu) - PCV`を求めます。
OSB適用後のコード/位相へIF係数で合成した`if_add_m`を一度加えます。
wind-upはさらに位相だけへ加えます。これは固定ローカル軸の受信機に対する一次の
視線投影モデルで、出力XYZはmarkerの座標です。位置微分や厳密なARP/APC座標での
送信時刻・Sagnac・大気の再評価は行いません。大きな設置オフセット向けの厳密モデルではありません。
`PppMeasurement.receiver_antenna`には製品名・型名・serial・元行・周波数、
PCO投影/PCV・各周波数とIFの加算値を残します。`observation`と生GF/MWは補正前です。

補正は各反復と最終後残差で再計算し、角度の滑らかな変化だけではアークを切りません。
校正はフィルターの寿命中に固定し、別の校正・設置条件への切替には新しいフィルターを
作成します。期限切れや周波数欠測はエポック更新を棄却し、数値平均/共分散を維持して
位相継続を切ります。重複/逆順エポックは状態を変えません。

実際の元C++ ANTEX受信機読込とCLAS受信機投影関数、通常PPPのNOAZI補間式による80ケースを
比較しています。通常C++のIF経路はIF-PCOで位置を移し、受信機PCVを適用しないため、
このRustのIF-PCV接続は拡張です。通常PPP全体との同一性を主張しません。
独立C++で9衛星×64エポックの観測へ補正を注入した2入力は全てPPP FLOAT、各568コード/位相を
採用し、継続アークも一致しました。合成雑音なし/ありの最大3次元誤差は約1.90 mm/0.695 m、
最終誤差は約0.067 mm/2.67 mmです。衛星COM→APCは下記で接続しています。衛星PCV、方位PCV、校正切替、
実配布ANTEX/実データの精度・速度・メモリ比較は今後の対象です。

## 衛星ANTEX IF-PCOとCOM→APC

```bash
cargo run --release --locked --offline --example ppp -- satellite
cargo run --release --locked --offline --example ppp -- satellite noisy
```

`satellite`はCOMを宣言した合成SP3と、衛星PCO・受信機PCO/PCV・wind-upを含む
合成観測を選び、対応する補正を有効にします。CSV末尾に`satellite_antenna_satellites`と
`satellite_antenna_max_offset_m`を加えます。既存の19列は衛星補正を無効にした場合に一致します。

`Antex::satellite_calibration(PRN, SVN, time)`で、PRN・正確なSVN・受信GPS時刻の期間内の
唯一の校正を選びます。期間の両端は含み、SVNや校正の自動切替は行いません。
`SatelliteAntennaModel::new(SatelliteAntennaBinding, calibrations)`へ選択済み校正を渡し、
`PppConfig.satellite_antenna = Some(model)`と`expected_reference_point = CentreOfMass`を
設定します。校正・追尾コードの組はフィルターの寿命中に固定します。
初期実装はGPS PRN 1〜32、GPS 1/2/5の同じ二周波追尾コードの組を全対象衛星に使います。
空のSVN・校正重複・必要周波数の欠測・退化した周波数の組を拒否します。

bindingにはANTEX製品名・追尾コード・`ClockSource`・時計基準名を指定します。
呼出側が確立した「この校正はこの精密時計の基準と互換」という宣言です。
時計製品ヘッダーから互換性を自動推定しません。PPP選択の追尾コード・時計基準名/出所と
実際のSP3/CLK選択を照合します。APC入力へのPCO設定や、PCOなしのCOM設定は拒否します。

ANTEXの衛星North/East/Upを機体x/y/zとして使い、受信機のようにNorth/Eastを交換しません。
IF係数で合成した機体PCOを、公称yawの`z = -COM/|COM|`、
`y = normalize(z × (Sun-COM))`、`x = y × z`でECEFへ回転して、送信COM位置へ一度加えます。
太陽位置はwind-upと同じ元の近似式による受信GPSラベルです。食/正午/深夜yawや
UT1/EOP、実際の姿勢・衛星PCVは適用しません。欠けた/退化した姿勢を前回値や別の軸で埋めません。

`PppModel.transmit`はCOMの位置・速度・公開時計・一度だけ補正した周期相対論を保持します。
`satellite_position_m`はAPC、`satellite_antenna`にはPRN/SVN・型名・元行・周波数/IF係数、
機体PCO・太陽/軸・ECEFオフセットとAPCを残します。APCで距離・標高角・Niell・設計行、
受信機補正とwind-upを評価し、PCOを観測へさらに加えません。
`PppSolution.reference_point`は入力SP3の宣言なのでCOMのままです。
元の適用順序に合わせ、COMで決めた送信時刻/速度/時計をPCO適用後に再計算しません。
PCOや姿勢の受信機座標微分は設計行に含めません。

期限切れ・必要な衛星校正の欠測・基準不一致・退化姿勢は、エポック更新を棄却して
数値平均/共分散を維持し、位相継続を切ります。衛星校正の失敗時に未補正COMを採用しません。
滑らかなPCO回転でアークを切らず、重複/逆順入力は既存状態を変えません。

元C++の実際の衛星ANTEX読込、IF係数と太陽ヘルパー、元PPPの回転式・精密幾何/Niellによる
81ケースで3軸・GPS L1/L2、逆順L2/L1、L1/L5を比較しました。
独立C++の完全な式から最後に一度丸めた観測は、雑音なし/あり各64エポックで全てPPP FLOAT、
各568コード/位相の採用とアーク継続を確認しました。合成最大/最終3次元誤差は
約1.78 mm/0.020 mm、雑音あり約0.696 m/2.61 mmです。
実配布ANTEXと時計基準の確認、実データ精度、元PPPProcessor全体の収束/状態/性能比較は未完了です。

## 初期固体地球潮汐をPPPへ接続する

```bash
cargo run --release --locked --offline --example ppp -- legacy-tide
cargo run --release --locked --offline --example ppp -- legacy-tide noisy
```

`legacy-tide`は近似固体地球潮汐・衛星IF-PCO・受信機PCO/PCV・wind-upを含む合成入力を
選び、全補正を有効にします。CSVへ`solid_tide_satellites`と`solid_tide_displacement_m`
（変位のノルム、m）を加えます。潮汐無効時は既存21列が一致し、追加2列はゼロです。

`PppConfig.receiver_reference = ReceiverPositionReference::NominalMarker`と
`solid_earth_tide = Some(SolidEarthTideModel::LegacyLoveApproximateSunMoon)`を組にします。
潮汐無効時の組は`InstantaneousMarker`/`None`で、矛盾する組は拒否します。
`PppSolution.receiver_reference`は推定・出力XYZの基準を示します。
`NominalMarker`は選択モデルの変位を加える前の座標であり、標準のtide-free/zero-tide
地球座標系への変換や恒久潮汐規約を宣言するものではありません。

`tides::solid_earth_tide(time, nominal_marker, model)`は太陽・月のECEF位置、個々と合計の
変位、`instantaneous_marker_m = nominal_marker_m + displacement_m`を返します。
PPPは各反復の基準座標へ変位を一度加え、瞬時座標で送信時刻・Sagnac距離・標高角・
Niell・設計行・受信機アンテナ・wind-upを評価します。出力は基準座標のままです。
`PppModel.receiver_position_m`と`solid_earth_tide`に評価座標と来歴を残します。
観測への別の潮汐距離加算は行わず、生のGF/MW診断・LLIも変えません。
不正な位置・天体・時刻をゼロ補正に置き換えず、エラーを返します。
失敗エポックでは数値平均・共分散・状態時刻を維持し、位相継続を切ります。

元C++の実際の近似月位置12件、固定Love数の太陽/月変位252件を比較しています。
変位は絶対2e-12 m、潮汐を含む精密幾何54件は位置/距離/予測コード1e-6 m以内です。
独立C++で全補正を注入した観測は、雑音なし/あり各64エポックがすべてPPP FLOAT、
各568コード/位相を採用し、最終3次元誤差は0.000011410 / 0.002606282 mです。
これは合成真値に対する回帰結果です。丸めなしC++観測の対照は最大0.000416336 m、
最終0.000000109 mで、位相を1e-4 cyclesへ丸めたRINEX入力の初期最大0.003571701 m
と区別します。各入力の許容差・SHA-256・再生成手順は[fixture記録](tests/fixtures/README.md)にあります。

元通常PPPの既定は**IERS2010 Step-1+2**です。このモデルは`use_iers_solid_tide=false`
で選ぶ代替の次数2固定Love数（h2=0.6078、l2=0.0847）に対応します。
近似太陽・月・GMSTは元と同じGPS尺度の日付で、UTC/UT1/EOPや高精度天体暦を使いません。
IERS2010の自動天体入力とPPP接続、恒久潮汐規約の変換、海洋荷重、極潮汐、大気荷重は未実装です。
元PPPProcessor全体・実データ精度・速度/メモリの同等性は未確認です。

## IERS2010固体地球潮汐の計算部品

```bash
cargo run --release --locked --offline --example iers_tide
```

この例は元SOFA/IERS経路で生成した地球固定の太陽・月座標とTT/UT引数を読み、150件の
変位と瞬時座標をCSVへ出力します。Rustで天体位置や時刻を生成する例、PPPの出力ではありません。

`iers2010::IersTideArguments::new(tt_centuries, ut_hour)`には、J2000からのTTユリウス世紀と
0以上24未満のUT時を渡します。元ラッパーはこの潮汐式のUTをUTCで近似します。
GPS週内秒やUTC MJDをそのままTT引数として渡すことはできません。
`EarthFixedSunMoon::new(sun_m, moon_m)`には、同一エポック・受信機と同一の地球固定座標系で
評価した地心太陽・月座標（m）を渡します。ICRS→ITRS回転・尺度・時刻の整合は呼出側の責任で、
ベクトルの大きさだけから座標系を検証できません。
`iers2010::solid_earth_tide(nominal_marker_m, bodies, arguments)`は基準座標へ一度加える変位と
瞬時座標、入力来歴、次数2/3・Step-1日周/半日周/緯度依存・Step-2日周/長周期の各項を返します。

元のGinan由来Dehant本体に従い、地心緯度、IERS半径6378136.6 m、緯度依存Love数、次数3、
非弾性補正、31日周・5長周期成分を使います。恒久潮汐を変換するStep-3は元と同じく省略します。
不正・非有限入力はエラーです。正確に地軸上の局は経度が定まらず、元の除算で非有限になるため
`AxialStation`を返します。極近傍の局はC++比較に含めています。

IERS公開基準例（許容差1 µm）、元の実際の本体・各項505件と、元SOFA/IERSラッパーの
地球固定天体入力150件（変位の絶対許容差2e-12 m）を照合しています。
日付境界の連続性、TT/UTによる周波数補正、不正位置・天体・時刻も検証します。
これは計算部品の照合で、自動天体暦・歳差章動・EOP製品読込や、IERS2010を使う
PPPの精度・状態・性能の検証ではありません。`PppConfig.solid_earth_tide`は引き続き
近似Love数モデルのみを受け付けます。再生成と来歴は[fixture記録](tests/fixtures/README.md)を参照してください。

## IERS入力の時刻尺度と地球回転

`earth_rotation::LeapSecondTable::historical(product_id, valid_until_utc_mjd)`は1980-01-06以降の
既知の正のうるう秒を2017-01-01まで保持します。有効終端UTC MJD（整数日、含まない）は必須で、
それまで追加の遷移がないことを呼出側が確認する契約です。最新の告示を自動取得する機能ではありません。
更新済み系列は`LeapSecondTable::new`へ明示的に渡し、重複・逆順・負のうるう秒・不正基準を拒否します。
従来の`CalendarDateTime`によるUTC読込政策と`SystemTime`のGPSラベル表現は維持します。

`RotationTimeScales::from_gps(gps, &leaps, eop)`はGPS→TAI（+19秒）、TT（TAI+32.184秒）、
UTC、UT1を評価し、うるう秒表の製品名と期限を残します。
`UtcEpoch`は実際の86400/86401秒で日内を割る準MJDで、うるう秒中の`seconds_of_day`は86400以上です。
均一尺度の`UniformMjd`は日とSI秒を分け、小数秒や週跨ぎを保持します。
`EarthOrientationParams`はUT1−UTC秒とxp/yp角秒を明示入力し、範囲外の時刻をゼロ補正で通しません。
`tide_arguments_utc_approx()`は元IERSラッパーと同じUT≈UTC規約で、TT世紀とUTC準日内時を作ります。

`earth_rotation_angle(ut1)`はIAU2000 ERA、`polar_motion_matrix(tt, eop)`は極運動と
TIO locator s′の−47 µas/世紀を評価します。
`compose_celestial_to_terrestrial(&times, celestial_to_intermediate)`は呼出側のCIO歳差章動行列と
ERA・極運動を合成します。行列の有限値・直交性・正の行列式は確認しますが、入力の座標系・
時刻・モデルを数値だけから確認できません。IAU2006/2000A歳差章動はまだ自動生成しません。

元SOFAの正しいGPS→TAI→UTC経路、時刻尺度・ERA・極運動・CIO合成を558ケースで照合しています。
18回のうるう秒、挿入中・日境界・小数秒・週境界、合成EOPを含みます。
元の一回だけのGPS→UTC MJD近似はうるう秒日に準MJDと異なり、比較対象の330行で差があり、
MJD差の秒換算で最大約1秒です。この差は記録し、SOFAのTAI→UTCを基準にします。
これは時刻・回転部品の検証です。EOP読込、天体暦・歳差章動の自動評価、IERS2010のPPP接続、
実測精度と元ソルバー全体の同等性能は未完了です。

## 移植元との意図的な違い

- 不正な入力や数値的な非収束は`Result`で返します。地球中心の測地座標は
  一意に決まらないためエラーです。
- 時刻の加算・構築は正負両方向の週跨ぎを正規化し、週の範囲超過を拒否します。
- 時刻の`PartialEq`は厳密な比較です。許容誤差は`difference_seconds`で明示します。
  元の近似等値と順序比較の不整合を持ち込みません。
- `SystemTime`変換は元のエポックオフセット規約を明示したGPS尺度です。
  小数秒を保持しますが、**GPSとUTCのうるう秒変換は行いません**。
  f64週内秒と`SystemTime`のナノ秒精度の範囲で丸めが生じます。
- WGS-84離心率は元の定数ヘッダーの扁平率から求めた値を使います。
  元の座標ヘッダーの丸め値との微小差を許容誤差付きで検証しています。
- 未知の衛星系やPRNゼロは有効な`SatelliteId`として構築できません。
  衛星系ごとのPRN割り当て判定は今後の読込層で扱います。
- 放送暦の`toes`は明示的な値です。BDTのtoeゼロを未設定センチネルとして
  GPS尺度のtoeで置換しません。暦の選択は絶対的な経過時間と健康状態を確認します。
- 欠測位置・速度・共分散は`None`です。測位結果は、衛星数が十分でも有限の
  位置を持たなければ有効になりません。ソルバー固有の診断項目は未移植です。
- 初期FLOATはゼロ状態も明示的に有効にし、SPD共分散をCholeskyで解き、
  Joseph形式で更新します。線形更新は元の全状態force-active更新と比較します。
  反復は同じpriorから再線形化し、一つの観測を複数回独立な情報として数えません。
  元RTKProcessorの全設定・反復手順とは同一ではありません。
- 初期FIXは上表の明示的な検証条件を使い、FLOATに整数を戻しません。
  元C++の固定位置計算・コード位相ゲートは比較しますが、共分散は独立したEigenの
  Schur補完を基準に検証します。元のFIX共分散出力や全体のFIX政策との比較ではありません。
- 精密プロダクトは時刻系・欠測・品質境界を確認し、既定で外挿を拒否します。
  元C++の有限なゼロ座標の採用、1点補間、広い外挿・欠測を飛ばした補間を持ち込みません。
  SP3/CLK時計を混ぜず、CLKが選ばれた衛星ではその系列を一貫して使います。
- 精密送信時刻経路は元の一次近似と再補間を明示的に区別し、両端点の品質・範囲を検査します。
  精密コードSPPは呼出側の信号/時計基準補正を必須にして、放送暦TGDへの暗黙の切替を行いません。
- 二周波モデルは正確な追尾コードと補正基準を明示し、縮退した周波数対をエラーにします。
  元の係数関数が返す単周波係数`(1, 0)`へ切り替えません。
- 位相アークは採用した更新だけを履歴に反映し、欠測・棄却・基準/補正変更で保守的に継続を切ります。
  GF/MWの通常PPP閾値は比較しますが、元ソルバーの全アーク政策の再現ではありません。

## 次の実装範囲

1. RTK/PPP比較用の出力・基準結果・データセット管理
2. PPPのIERS2010潮汐に必要な天体暦・UTC/UT1/EOPを整え、海洋荷重・極潮汐、衛星PCV・姿勢・実プロダクトの補正基準対応を進める
3. GLONASSとマルチGNSS、補正・異常検知・FIX判定の移植
4. 元C++との実データ・状態遷移・速度比較で[最終目標](GOALS.md)を検証

コアはMITライセンス。移植元の著作権表示を[LICENSE](LICENSE)に保持しています。
LAMBDAのRTKLIB由来部分にはBSD 2-clauseの条件が適用され、元の表記を
[LICENSE-RTKLIB](LICENSE-RTKLIB)に保持しています。
Ginan由来のIERS2010 Rust移植にはApache 2.0が適用され、
[LICENSE-IERS2010](LICENSE-IERS2010)と[NOTICE-IERS2010](NOTICE-IERS2010)に条件・来歴・変更点を保持しています。
パッケージの表記は`MIT AND BSD-2-Clause AND Apache-2.0`です。
SOFAソースはRustクレートに含めず、C++基準生成時に移植元の未変更コピーを利用します。
