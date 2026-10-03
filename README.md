# JMA GRIB2 to NetCDF (QGIS plugin)

気象庁のランレングス圧縮GRIB2（データ表現テンプレート5.200）を、時刻次元を持つCF準拠のNetCDF-4に変換するQGISプラグインです。
追加ライブラリは不要です（numpy とQGIS同梱のGDALだけで動きます）。wgrib2・ecCodesも使いません。

## 背景

- GDALのGRIBドライバは、テンプレート5.200に対応していません（最新版でも未対応）。
- ecCodesは5.200に対応していますが、解析雨量で使われる気象庁独自のプロダクト定義テンプレート4.50008が定義に入っていないため、読み込みが途中で止まります。
- そこでこのプラグインは、必要なセクションを自前で解析し、ランレングス展開をnumpyでベクトル化して処理します。全国1km格子（約860万点）1時刻あたり約0.1秒です。

## 対応データ（v0.1）

| 項目 | 対応 |
|---|---|
| データ表現 | 5.200（ランレングス圧縮）のみ。それ以外のファイルはスキップしてログに残す |
| 格子 | 3.0（等緯度経度） |
| 時刻 | テンプレート4.8／4.50008は「積算期間の終わり」、それ以外は「参照時刻＋予報時間」。すべてUTC |
| 検証済み | 解析雨量（`*_SRF_GPV_Ggis1km_Prr60lv_ANAL_grib2.bin`）。ecCodesのデコード結果と全点一致 |

## 使い方

- メニューの **ラスタ → JMA GRIB2 to NetCDF**、またはProcessingツールボックスの「気象庁GRIB2（ランレングス）→ NetCDF」から実行します。
- 入力はフォルダ単位です（ファイル名パターン指定、サブフォルダの検索に対応）。期間を指定する場合はUTCで入力します。
- 出力は `precip(time, lat, lon)`（float32、欠測は -9999）で、`time_bnds`（積算期間）と、測地系JGD2000の座標系情報（EPSG:4612）を含みます。
- QGISでは「1バンド＝1時刻」のラスタとして開きます。

## 開発

```
core/        QGIS非依存（grib2.py: デコーダ, converter.py: NetCDF書き出し）
algorithms/  Processingアルゴリズム
tests/       pytest（合成データでの往復テスト。実データは JMA_GRIB2_SAMPLE=<file> で追加テスト）
```

```
python -m pytest jma_grib2nc/tests
flake8 --max-line-length=100 jma_grib2nc
bandit -r jma_grib2nc -x jma_grib2nc/tests
```

License: GPL-2.0-or-later
