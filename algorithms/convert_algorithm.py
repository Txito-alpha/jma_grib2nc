# -*- coding: utf-8 -*-
"""Processing algorithm: JMA GRIB2 (template 5.200) -> NetCDF."""
import os

from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingOutputRasterLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDateTime,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterString,
    QgsRasterLayer,
)
from qgis.PyQt.QtCore import QCoreApplication

from ..core.converter import ConversionError, Reporter, convert, find_files


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
            "・欠測は -9999（NoData）\n\n"
            "Converts JMA run-length packed GRIB2 files (template 5.200) in a folder "
            "into a single CF-compliant NetCDF-4 file with a time dimension. "
            "Times are UTC."
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

        files = find_files(folder, pattern, recursive)
        feedback.pushInfo(f"input folder: {folder} (pattern {pattern}, "
                          f"recursive={recursive}) -> {len(files)} file(s)")
        if start or end:
            feedback.pushInfo(f"period (UTC): {start or '-'} .. {end or '-'}")
        if not files:
            raise QgsProcessingException(self.tr("パターンに一致するファイルがありません"))

        try:
            result = convert(files, out_path, _FeedbackReporter(feedback),
                             start=start, end=end, zlevel=zlevel)
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

    def groupId(self):
        return ""

    def group(self):
        return ""

    def helpUrl(self):
        return "https://github.com/Txito-alpha/jma_grib2nc"
