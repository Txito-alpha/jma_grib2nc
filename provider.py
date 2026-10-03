# -*- coding: utf-8 -*-
import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .algorithms.convert_algorithm import ConvertJmaGrib2Algorithm


class JmaGrib2NcProvider(QgsProcessingProvider):
    def loadAlgorithms(self):
        self.addAlgorithm(ConvertJmaGrib2Algorithm())

    def id(self):
        return "jmagrib2nc"

    def name(self):
        return "JMA GRIB2 to NetCDF"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))
