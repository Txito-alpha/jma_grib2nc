# -*- coding: utf-8 -*-
"""Default symbology applied to layers added by the plugin."""
from qgis.core import (
    QgsColorRampShader,
    QgsGradientColorRamp,
    QgsMeshDatasetIndex,
    QgsRectangle,
    QgsStyle,
)
from qgis.PyQt.QtGui import QColor

# Mesh "Contours" (scalar) style for precipitation
PRECIP_DATASET_NAME = "60-minute accumulated precipitation"
PRECIP_MIN = 0.0
PRECIP_MAX = 80.0
PRECIP_CLASSES = 17          # equal interval -> 0, 5, 10, ... 80 mm
PRECIP_RAMP = "Spectral"     # inverted: blue (light rain) -> red (heavy rain)


def _precip_ramp():
    ramp = QgsStyle.defaultStyle().colorRamp(PRECIP_RAMP)
    if ramp is None:  # style database without the ramp: same end colours as Spectral
        ramp = QgsGradientColorRamp(QColor("#d7191c"), QColor("#2b83ba"))
    ramp.invert()
    return ramp


def precip_color_ramp_shader():
    """Interpolated colour ramp 0-80 mm, equal interval 17 classes, value 0 transparent."""
    shader = QgsColorRampShader(PRECIP_MIN, PRECIP_MAX, _precip_ramp(),
                                QgsColorRampShader.Interpolated,
                                QgsColorRampShader.EqualInterval)
    shader.classifyColorRamp(PRECIP_CLASSES, -1, QgsRectangle(), None)
    items = shader.colorRampItemList()
    for item in items:
        if item.value == 0:
            c = QColor(item.color)
            c.setAlpha(0)
            item.color = c
    shader.setColorRampItemList(items)
    return shader


def find_dataset_group(layer, name):
    """Index of the mesh dataset group with this name (first group if not found)."""
    groups = layer.datasetGroupsIndexes()
    for g in groups:
        if layer.datasetGroupMetadata(QgsMeshDatasetIndex(g)).name() == name:
            return g
    return groups[0] if groups else -1


def apply_mesh_precip_style(layer):
    """Show the precipitation dataset as contours (scalar) with the plugin's default style.

    Returns the dataset group index that was styled, or -1 when the layer has none.
    """
    group = find_dataset_group(layer, PRECIP_DATASET_NAME)
    if group < 0:
        return -1
    rs = layer.rendererSettings()
    sc = rs.scalarSettings(group)
    sc.setClassificationMinimumMaximum(PRECIP_MIN, PRECIP_MAX)
    sc.setColorRampShader(precip_color_ramp_shader())
    rs.setScalarSettings(group, sc)
    rs.setActiveScalarDatasetGroup(group)
    layer.setRendererSettings(rs)
    layer.triggerRepaint()
    return group
