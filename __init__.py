# -*- coding: utf-8 -*-
"""JMA GRIB2 to NetCDF - QGIS plugin entry point."""
import os
import sys

# On Windows QGIS, sys.stderr/stdout can be None, which turns real import errors
# into misleading AttributeErrors. Restore harmless streams before anything else.
if sys.stderr is None:
    sys.stderr = open(os.devnull, "w")
if sys.stdout is None:
    sys.stdout = open(os.devnull, "w")


def classFactory(iface):  # noqa: N802 (QGIS API name)
    from .plugin import JmaGrib2NcPlugin
    return JmaGrib2NcPlugin(iface)
