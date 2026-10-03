# -*- coding: utf-8 -*-
"""Headless tests for the core package (no QGIS needed).

Run from the folder that contains jma_grib2nc:  python -m pytest jma_grib2nc/tests
Set JMA_GRIB2_SAMPLE to a real JMA GRIB2 file to also run the real-data test.
"""
import os
import struct
from datetime import datetime

import numpy as np
import pytest

from jma_grib2nc.core import grib2
from jma_grib2nc.core.converter import (ConversionError, Reporter, collect, convert,
                                        subset_window)


def _s32(v):
    return struct.pack(">I", (abs(v) | 0x80000000) if v < 0 else v)


def encode_rle(levels_idx, nbit=8, vmax=None):
    """Reference encoder for template 5.200 (used to test the decoder)."""
    vmax = vmax if vmax is not None else int(levels_idx.max())
    lngu = (1 << nbit) - 1 - vmax
    out, i, n = [], 0, len(levels_idx)
    while i < n:
        lv = int(levels_idx[i])
        j = i + 1
        while j < n and levels_idx[j] == lv:
            j += 1
        out.append(lv)
        r = j - i - 1
        while r > 0:
            out.append(vmax + 1 + r % lngu)
            r //= lngu
        i = j
    if nbit == 8:
        return bytes(out), vmax
    bits = "".join(format(v, f"0{nbit}b") for v in out)
    bits += "0" * (-len(bits) % 8)
    return bytes(int(bits[k:k + 8], 2) for k in range(0, len(bits), 8)), vmax


def build_message(lv_grid, level_values, dscale=1, nbit=8,
                  t=datetime(2016, 8, 30, 9), interval=60, scan=0):
    nj, ni = lv_grid.shape
    data, vmax = encode_rle(lv_grid.ravel(), nbit)
    s1 = struct.pack(">IBHHBBBHBBBBBBB", 21, 1, 34, 0, 2, 1, 0, t.year, t.month,
                     t.day, t.hour, t.minute, 0, 0, 2)
    s3 = (struct.pack(">IBBIBBH", 72, 3, 0, ni * nj, 0, 0, 0) + bytes([4]) + b"\xff" * 15
          + struct.pack(">II", ni, nj) + b"\x00" * 8
          + _s32(48_000_000) + _s32(140_000_000) + b"\x30"
          + _s32(48_000_000 - (nj - 1) * 10000) + _s32(140_000_000 + (ni - 1) * 10000)
          + struct.pack(">II", 10000, 10000) + bytes([scan]))
    assert len(s3) == 72
    body4 = (bytes([1, 200, 0, 150, 255]) + b"\x00\x00\x0a" + bytes([0]) + _s32(-interval)
             + bytes([1, 255]) + b"\xff" * 4 + bytes([255, 255]) + b"\xff" * 4
             + struct.pack(">HBBBBB", t.year, t.month, t.day, t.hour, t.minute, 0)
             + bytes([1]) + b"\x00" * 4 + bytes([1, 2, 0]) + struct.pack(">I", interval)
             + bytes([0]) + b"\x00" * 4)
    s4 = struct.pack(">IBHH", 9 + len(body4), 4, 0, 50008) + body4
    lv_bytes = b"".join(struct.pack(">H", v) for v in level_values)
    s5 = (struct.pack(">IBIHB", 17 + len(lv_bytes), 5, ni * nj, 200, nbit)
          + struct.pack(">HHB", vmax, len(level_values), dscale) + lv_bytes)
    s6 = struct.pack(">IBB", 6, 6, 255)
    s7 = struct.pack(">IB", 5 + len(data), 7) + data
    body = s1 + s3 + s4 + s5 + s6 + s7 + b"7777"
    return b"GRIB\x00\x00\x00\x02" + struct.pack(">Q", 16 + len(body)) + body


@pytest.fixture
def sample(tmp_path):
    rng = np.random.default_rng(0)
    nj, ni = 40, 30
    lv = np.zeros((nj, ni), dtype=np.int64)          # 0 = missing
    lv[5:35, 3:27] = 1                                # 0 mm
    lv[10:20, 5:15] = rng.integers(1, 6, size=(10, 10))
    lv[25, 20] = 5
    level_values = [0, 4, 10, 20, 500]                # /10 -> 0, 0.4, 1, 2, 50 mm
    expect = np.concatenate([[np.nan], np.array(level_values) / 10.0])[lv]
    files = []
    for k, hour in enumerate((9, 10)):
        p = tmp_path / f"msg{k}.bin"
        p.write_bytes(build_message(lv, level_values, t=datetime(2016, 8, 30, hour)))
        files.append(str(p))
    return files, expect


def test_decode_roundtrip(sample):
    files, expect = sample
    m = grib2.scan_file(files[0])[0]
    assert m.drt == 200 and m.pdt == 50008
    assert m.valid_time == datetime(2016, 8, 30, 9)
    assert m.param_key == (0, 1, 200, 1, 60)
    got = grib2.decode(m)
    np.testing.assert_array_equal(np.isnan(got), np.isnan(expect))
    np.testing.assert_allclose(got[~np.isnan(got)], expect[~np.isnan(expect)])


def test_long_runs_and_nbit():
    lv = np.array([1] * 100000 + [2] + [1] * 3, dtype=np.int64)
    for nbit in (4, 8, 12):
        data, vmax = encode_rle(lv, nbit=nbit, vmax=2)
        s5 = (struct.pack(">IBIHB", 21, 5, lv.size, 200, nbit)
              + struct.pack(">HHB", vmax, 2, 0) + struct.pack(">HH", 7, 9))
        s7 = struct.pack(">IB", 5 + len(data), 7) + data
        out = grib2.decode_run_length(s5, s7, lv.size)
        np.testing.assert_array_equal(out, np.where(lv == 1, 7.0, 9.0))


def test_grid_coordinates(sample):
    files, _ = sample
    g = grib2.scan_file(files[0])[0].grid
    assert g.lats()[0] == pytest.approx(48.0) and g.lats()[-1] == pytest.approx(47.61)
    assert g.lons()[-1] == pytest.approx(140.29)


def test_collect_skips_bad_files(sample, tmp_path):
    files, _ = sample
    bad = tmp_path / "bad.bin"
    bad.write_bytes(b"not grib")
    warnings = []

    class R(Reporter):
        def info(self, msg):
            pass

        def warn(self, msg):
            warnings.append(msg)

    grid, variables, skipped = collect(files + [str(bad)], R())
    assert skipped == 1 and any("bad.bin" in w for w in warnings)
    assert list(variables.values())[0].name == "precip"


def test_convert_netcdf(sample, tmp_path):
    gdal = pytest.importorskip("osgeo.gdal")
    files, expect = sample
    out = tmp_path / "出力" / "テスト.nc"
    res = convert(files, str(out), start=datetime(2016, 8, 30, 10))
    assert res["times"] == 1 and res["variables"] == ["precip"]
    ds = gdal.Open(f'NETCDF:"{out}":precip')
    a = ds.GetRasterBand(1).ReadAsArray()
    np.testing.assert_allclose(np.where(np.isnan(expect), -9999, expect), a, rtol=1e-6)
    gt = ds.GetGeoTransform()
    assert gt[0] == pytest.approx(139.995) and gt[3] == pytest.approx(48.005)


def test_subset_window(sample):
    files, _ = sample
    g = grib2.scan_file(files[0])[0].grid   # lon 140.00..140.29, lat 48.00..47.61, 0.01 deg
    # box inside the grid: every cell touching the box is taken
    assert subset_window(g, (140.05, 47.80, 140.10, 47.90)) == (10, 21, 5, 11, False)
    # box edges exactly on cell edges: neighbouring cells are not added
    assert subset_window(g, (140.045, 47.795, 140.105, 47.905)) == (10, 21, 5, 11, False)
    # partly outside -> clipped to the grid
    assert subset_window(g, (139.0, 47.0, 140.02, 47.62)) == (38, 40, 0, 3, True)
    with pytest.raises(ConversionError):
        subset_window(g, (150.0, 30.0, 151.0, 31.0))
    with pytest.raises(ConversionError):
        subset_window(g, (140.2, 47.7, 140.1, 47.8))


def test_subset_window_rounded_corners():
    # JMA 1 km grid: corner coordinates rounded to 1e-6 deg (1/120 deg cells)
    g = grib2.Grid(ni=2560, nj=3360, lat1=47.995833, lon1=118.00625, lat2=20.004167,
                   lon2=149.99375, dx=0.0125, dy=0.008333, scan=0, earth_shape=4)
    # 45.6N / 41.3N / 139.3E / 145.9E are exact cell edges: no extra row or column
    assert subset_window(g, (139.3, 41.3, 145.9, 45.6)) == (288, 804, 1704, 2232, False)
    # a point: the cell containing it (both neighbours when on a cell edge), widened to at
    # least 2 x 2 cells because GDAL cannot georeference 1-cell-wide netCDF rasters
    assert subset_window(g, (141.355, 43.06, 141.355, 43.06)) == (592, 594, 1868, 1870, False)
    assert subset_window(g, (141.35, 43.06, 141.35, 43.06)) == (592, 594, 1867, 1869, False)
    # at the last row / column the neighbour is taken on the inner side
    assert subset_window(g, (149.999, 20.001, 149.999, 20.001)) == (3358, 3360, 2558, 2560, False)


def test_convert_subset(sample, tmp_path):
    gdal = pytest.importorskip("osgeo.gdal")
    files, expect = sample
    out = tmp_path / "subset.nc"
    res = convert(files, str(out), bbox=(140.05, 47.80, 140.10, 47.90))
    assert res["window"] == (10, 21, 5, 11)
    ds = gdal.Open(f'NETCDF:"{out}":precip')
    assert (ds.RasterXSize, ds.RasterYSize, ds.RasterCount) == (6, 11, 2)
    a = ds.GetRasterBand(1).ReadAsArray()
    sub = expect[10:21, 5:11]
    np.testing.assert_allclose(np.where(np.isnan(sub), -9999, sub), a, rtol=1e-6)
    gt = ds.GetGeoTransform()
    assert gt[0] == pytest.approx(140.045) and gt[3] == pytest.approx(47.905)
    assert gt[1] == pytest.approx(0.01) and gt[5] == pytest.approx(-0.01)


@pytest.mark.skipif(not os.environ.get("JMA_GRIB2_SAMPLE"), reason="JMA_GRIB2_SAMPLE not set")
def test_real_file():
    m = grib2.scan_file(os.environ["JMA_GRIB2_SAMPLE"])[0]
    d = grib2.decode(m)
    assert d.shape == (m.grid.nj, m.grid.ni)
    assert np.nanmin(d) >= 0
