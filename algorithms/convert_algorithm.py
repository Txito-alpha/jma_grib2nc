# -*- coding: utf-8 -*-
"""Processing algorithm: JMA GRIB2 (template 5.200) -> NetCDF."""
import os

from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeatureRequest,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingOutputRasterLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDateTime,
    QgsProcessingParameterExtent,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterString,
    QgsRasterLayer,
    QgsRectangle,
)
from qgis.PyQt.QtCore import QCoreApplication

from ..core.converter import ConversionError, Reporter, convert, find_files

# JMA GRIB2 grids use the GRS80 ellipsoid (JGD2000). The difference to WGS 84 / JGD2011
# is far below the 1 km cell size, so extents are transformed to EPSG:4612.
GRID_CRS = "EPSG:4612"


class _FeedbackReporter(Reporter):
    def __init__(self, feedback):
        self.fb = feedback

    def info(self, msg):
        self.fb.pushInfo(msg)

    def warn(self, msg):
        self.fb.pushWarning(msg)

    def progress(self, percent):
        self.fb.setProgress(percent)

    def canceled(self):
        return self.fb.isCanceled()


def _to_datetime(qdt):
    if qdt is None or not qdt.isValid():
        return None
    return qdt.toPyDateTime().replace(tzinfo=None)


class ConvertJmaGrib2Algorithm(QgsProcessingAlgorithm):
    INPUT_FOLDER = "INPUT_FOLDER"
    PATTERN = "PATTERN"
    RECURSIVE = "RECURSIVE"
    START = "START"
    END = "END"
    EXTENT = "EXTENT"
    EXTENT_FEATURES = "EXTENT_FEATURES"
    ZLEVEL = "ZLEVEL"
    LOAD = "LOAD"
    OUTPUT = "OUTPUT"
    OUTPUT_LAYER = "OUTPUT_LAYER"

    def tr(self, text):
        return QCoreApplication.translate("ConvertJmaGrib2Algorithm", text)

    def createInstance(self):
        return ConvertJmaGrib2Algorithm()

    def name(self):
        return "convert_jma_grib2_rle"

    def displayName(self):
        return self.tr("気象庁GRIB2（ランレングス）→ NetCDF")

    def shortHelpString(self):
        return self.tr(
            "気象庁のランレングス圧縮GRIB2（データ表現テンプレート5.200：解析雨量、"
            "降水短時間予報、降水ナウキャスト等）をフォルダ単位で読み込み、"
            "時刻次元を持つCF準拠のNetCDF-4に変換します。\n\n"
            "・追加ライブラリ不要（numpy と QGIS同梱のGDALのみ）\n"
            "・時刻はUTC（積算値は積算期間の終わりの時刻）。期間指定もUTCで入力\n"
            "・読めないファイル、未対応形式、格子の異なるファイルはスキップしてログに記録\n"
            "・欠測は -9999（NoData）\n"
            "・出力範囲：範囲（地図キャンバス／レイヤの範囲／描画など）または地物"
            "（「選択地物のみ」可）の範囲で切り出し可能。範囲にかかる格子セルをすべて含むよう"
            "セル境界に合わせて外側に広げ、座標は元の格子のまま。CRSは自動変換。"
            "両方空欄なら全域\n\n"
            "Converts JMA run-length packed GRIB2 files (template 5.200) in a folder "
            "into a single CF-compliant NetCDF-4 file with a time dimension. "
            "Times are UTC. Optionally limit the output to an extent or to the extent of "
            "(selected) features; the window is snapped outward to whole grid cells."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            self.INPUT_FOLDER, self.tr("入力フォルダ（GRIB2）"),
            behavior=QgsProcessingParameterFile.Folder))
        self.addParameter(QgsProcessingParameterString(
            self.PATTERN, self.tr("ファイル名パターン"), defaultValue="*.bin"))
        self.addParameter(QgsProcessingParameterBoolean(
            self.RECURSIVE, self.tr("サブフォルダも検索する"), defaultValue=False))
        self.addParameter(QgsProcessingParameterDateTime(
            self.START, self.tr("開始日時（UTC、空欄で制限なし）"), optional=True))
        self.addParameter(QgsProcessingParameterDateTime(
            self.END, self.tr("終了日時（UTC、空欄で制限なし）"), optional=True))
        self.addParameter(QgsProcessingParameterExtent(
            self.EXTENT, self.tr("出力範囲（空欄で全域）"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.EXTENT_FEATURES, self.tr("地物の範囲で切り出す（空欄で使わない）"),
            [QgsProcessing.TypeVectorAnyGeometry], optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.ZLEVEL, self.tr("圧縮レベル（0=無圧縮, 1-9）"),
            QgsProcessingParameterNumber.Integer, defaultValue=4, minValue=0, maxValue=9))
        self.addParameter(QgsProcessingParameterBoolean(
            self.LOAD, self.tr("変換後に地図に追加する"), defaultValue=True))
        self.addParameter(QgsProcessingParameterFileDestination(
            self.OUTPUT, self.tr("出力NetCDF"), fileFilter="NetCDF (*.nc)"))
        self.addOutput(QgsProcessingOutputRasterLayer(
            self.OUTPUT_LAYER, self.tr("変換結果レイヤ")))

    def processAlgorithm(self, parameters, context, feedback):
        folder = self.parameterAsFile(parameters, self.INPUT_FOLDER, context)
        pattern = self.parameterAsString(parameters, self.PATTERN, context) or "*"
        recursive = self.parameterAsBoolean(parameters, self.RECURSIVE, context)
        start = _to_datetime(self.parameterAsDateTime(parameters, self.START, context)) \
            if parameters.get(self.START) else None
        end = _to_datetime(self.parameterAsDateTime(parameters, self.END, context)) \
            if parameters.get(self.END) else None
        zlevel = self.parameterAsInt(parameters, self.ZLEVEL, context)
        load = self.parameterAsBoolean(parameters, self.LOAD, context)
        out_path = self.parameterAsFileOutput(parameters, self.OUTPUT, context)
        if not out_path.lower().endswith(".nc"):
            out_path += ".nc"

        if not folder or not os.path.isdir(folder):
            raise QgsProcessingException(self.tr("入力フォルダが見つかりません: ") + str(folder))
        if start and end and start > end:
            raise QgsProcessingException(self.tr("開始日時が終了日時より後になっています"))

        bbox = self._bbox(parameters, context, feedback)

        files = find_files(folder, pattern, recursive)
        feedback.pushInfo(f"input folder: {folder} (pattern {pattern}, "
                          f"recursive={recursive}) -> {len(files)} file(s)")
        if start or end:
            feedback.pushInfo(f"period (UTC): {start or '-'} .. {end or '-'}")
        if not files:
            raise QgsProcessingException(self.tr("パターンに一致するファイルがありません"))

        try:
            result = convert(files, out_path, _FeedbackReporter(feedback),
                             start=start, end=end, zlevel=zlevel, bbox=bbox)
        except ConversionError as e:
            raise QgsProcessingException(str(e))
        if result.get("canceled"):
            feedback.pushWarning("canceled; partial output removed")
            return {}

        outputs = {self.OUTPUT: out_path}
        if result["skipped"]:
            feedback.pushWarning(f"{result['skipped']} file/message(s) skipped - see log above")
        var = result["variables"][0]
        uri = f'NETCDF:"{out_path}":{var}'
        outputs[self.OUTPUT_LAYER] = uri
        if load:
            name = f"{os.path.splitext(os.path.basename(out_path))[0]} ({var})"
            layer = QgsRasterLayer(uri, name, "gdal")
            if layer.isValid():
                context.temporaryLayerStore().addMapLayer(layer)
                context.addLayerToLoadOnCompletion(
                    layer.id(),
                    QgsProcessingContext.LayerDetails(name, context.project(), self.OUTPUT_LAYER))
                feedback.pushInfo(f"layer added: {name}, {layer.bandCount()} band(s) "
                                  "(1 band = 1 time step)")
            else:
                feedback.pushWarning("output written but could not be opened as a layer: " + uri)
        return outputs

    def _bbox(self, parameters, context, feedback):
        """Requested output extent as (west, south, east, north) in GRID_CRS, or None."""
        grid_crs = QgsCoordinateReferenceSystem(GRID_CRS)
        use_extent = bool(parameters.get(self.EXTENT))
        use_features = bool(parameters.get(self.EXTENT_FEATURES))
        if use_extent and use_features:
            raise QgsProcessingException(
                self.tr("「出力範囲」と「地物の範囲」はどちらか一方だけ指定してください"))
        if use_features:
            source = self.parameterAsSource(parameters, self.EXTENT_FEATURES, context)
            if source is None:
                raise QgsProcessingException(self.tr("範囲に使う地物レイヤを開けません"))
            rect = None
            n = 0
            request = QgsFeatureRequest().setNoAttributes()
            for f in source.getFeatures(request):
                if f.hasGeometry() and not f.geometry().isEmpty():
                    bb = f.geometry().boundingBox()
                    if rect is None:
                        rect = QgsRectangle(bb)
                    else:
                        rect.combineExtentWith(bb)
                    n += 1
            if rect is None:
                raise QgsProcessingException(self.tr("範囲に使う地物（ジオメトリ）がありません"))
            src_crs = source.sourceCrs()
            if not src_crs.isValid():
                feedback.pushWarning(f"feature layer has no CRS; assumed {GRID_CRS}")
                src_crs = grid_crs
            if src_crs != grid_crs:
                tr = QgsCoordinateTransform(src_crs, grid_crs, context.transformContext())
                rect = tr.transformBoundingBox(rect)
            feedback.pushInfo(f"extent from {n} feature(s) ({src_crs.authid()})")
        elif use_extent:
            rect = self.parameterAsExtent(parameters, self.EXTENT, context, grid_crs)
            src_crs = self.parameterAsExtentCrs(parameters, self.EXTENT, context)
            feedback.pushInfo(f"extent given in {src_crs.authid() or 'unknown CRS'}")
        else:
            return None
        if rect.isNull():
            raise QgsProcessingException(self.tr("出力範囲が空です"))
        feedback.pushInfo(f"requested extent ({GRID_CRS}): "
                          f"lon {rect.xMinimum():.5f}..{rect.xMaximum():.5f}, "
                          f"lat {rect.yMinimum():.5f}..{rect.yMaximum():.5f}")
        return (rect.xMinimum(), rect.yMinimum(), rect.xMaximum(), rect.yMaximum())

    def groupId(self):
        return ""

    def group(self):
        return ""

    def helpUrl(self):
        return "https://github.com/Txito-alpha/jma_grib2nc"
