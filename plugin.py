# -*- coding: utf-8 -*-
import os

from qgis.core import QgsApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction

from .provider import JmaGrib2NcProvider

ALGORITHM_ID = "jmagrib2nc:convert_jma_grib2_rle"
MENU = "JMA GRIB2 to NetCDF"


class JmaGrib2NcPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.provider = None
        self.action = None

    def initProcessing(self):
        self.provider = JmaGrib2NcProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()
        icon = QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))
        self.action = QAction(icon, "気象庁GRIB2（ランレングス）→ NetCDF…",
                              self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu(MENU, self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        if self.action:
            self.iface.removePluginMenu(MENU, self.action)
            self.iface.removeToolBarIcon(self.action)
            self.action = None
        if self.provider:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None

    def run(self):
        import processing
        processing.execAlgorithmDialog(ALGORITHM_ID, {})
