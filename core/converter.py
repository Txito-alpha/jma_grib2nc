# -*- coding: utf-8 -*-
"""Convert JMA GRIB2 (template 5.200) files to a CF-compliant NetCDF-4 file.

QGIS independent: depends on numpy and GDAL (netCDF multidimensional API, GDAL >= 3.1).
Problems with individual files are logged and the file is skipped; the run only
fails when nothing usable is left.
"""
import os
import shutil
import tempfile
from collections import OrderedDict
from datetime import datetime, timezone
from typing import Dict, List, Optional, Sequence

import numpy as np
from osgeo import gdal

from .grib2 import Grib2Error, MessageInfo, decode, scan_file

FILL_VALUE = -9999.0
EPOCH = datetime(1970, 1, 1)
TOOL_NAME = "JMA GRIB2 to NetCDF (QGIS plugin)"

_EARTH = {  # code table 3.2 -> (semi-major axis, inverse flattening, CRS name, EPSG)
    4: (6378137.0, 298.257222101, "JGD2000", 4612),
    5: (6378137.0, 298.257223563, "WGS 84", 4326),
}


class Reporter:
    """Logging / progress hooks. The Processing algorithm overrides these."""

    def info(self, msg: str) -> None:
        print(msg)

    def warn(self, msg: str) -> None:
        print("WARNING: " + msg)

    def progress(self, percent: float) -> None:
        pass

    def canceled(self) -> bool:
        return False


class ConversionError(Exception):
    pass


class Variable:
    def __init__(self, key, name, attrs):
        self.key = key
        self.name = name
        self.attrs = attrs
        self.messages: Dict[datetime, MessageInfo] = {}


def _describe(key) -> dict:
    discipline, category, number, stat, interval = key
    if discipline == 0 and category == 1:
        if stat == 1 and interval:
            return dict(base="precip", units="mm", standard_name="precipitation_amount",
                        long_name=f"{interval}-minute accumulated precipitation",
                        cell_methods=f"time: sum (interval: {interval} minutes)")
        return dict(base="precip_rate", units="mm h-1",
                    long_name="precipitation intensity")
    return dict(base=f"var_d{discipline}_c{category}_n{number}", units="1",
                long_name=f"JMA GRIB2 parameter discipline={discipline} "
                          f"category={category} number={number}")


def find_files(folder: str, pattern: str = "*.bin", recursive: bool = False) -> List[str]:
    import fnmatch
    hits = []
    for root, dirs, files in os.walk(folder):
        dirs.sort()
        for name in sorted(files):
            if fnmatch.fnmatch(name.lower(), pattern.lower()):
                hits.append(os.path.join(root, name))
        if not recursive:
            break
    return hits


def collect(paths: Sequence[str], reporter: Reporter,
            start: Optional[datetime] = None, end: Optional[datetime] = None):
    """Scan headers, filter by valid time (UTC) and group messages into variables."""
    grid = None
    variables: "OrderedDict[tuple, Variable]" = OrderedDict()
    n_skipped = 0
    for i, path in enumerate(paths):
        if reporter.canceled():
            return None, variables, n_skipped
        reporter.progress(10.0 * (i + 1) / max(len(paths), 1))
        name = os.path.basename(path)
        try:
            msgs = scan_file(path)
        except (Grib2Error, OSError, ValueError) as e:
            reporter.warn(f"skip {name}: {e}")
            n_skipped += 1
            continue
        for m in msgs:
            if "warning" in m.extra:
                reporter.warn(f"{name}: {m.extra['warning']}")
            if m.drt != 200:
                reporter.warn(f"skip {name} (message @{m.offset}): data template 5.{m.drt} "
                              "is not supported (only 5.200 run-length products)")
                n_skipped += 1
                continue
            if grid is None:
                grid = m.grid
                reporter.info(f"grid: {grid.ni} x {grid.nj}, "
                              f"lat {grid.lat2}..{grid.lat1}, lon {grid.lon1}..{grid.lon2}")
            elif m.grid != grid:
                reporter.warn(f"skip {name}: grid differs from the first file")
                n_skipped += 1
                continue
            if (start and m.valid_time < start) or (end and m.valid_time > end):
                continue
            var = variables.get(m.param_key)
            if var is None:
                d = _describe(m.param_key)
                var = Variable(m.param_key, d.pop("base"), d)
                variables[m.param_key] = var
            if m.valid_time in var.messages:
                reporter.warn(f"duplicate {var.name} at {m.valid_time:%Y-%m-%d %H:%M} UTC "
                              f"in {name}; first one kept")
                continue
            var.messages[m.valid_time] = m
    # unique variable names
    seen: Dict[str, int] = {}
    for var in variables.values():
        if var.name in seen:
            seen[var.name] += 1
            var.name = f"{var.name}_{seen[var.name]}"
        else:
            seen[var.name] = 0
    return grid, variables, n_skipped


def subset_window(grid, bbox, tol: Optional[float] = None):
    """Grid index window covering a lon/lat bounding box.

    bbox = (west, south, east, north) in the grid's geographic CRS (degrees).
    Every cell whose area intersects the box is selected, so the output always
    covers the requested extent and keeps the original cell centres.
    Returns (j0, j1, i0, i1, clipped) with half-open ranges in output order
    (rows north -> south); clipped is True when the box extends beyond the grid.
    Raises ConversionError when the box does not overlap the grid.

    tol defaults to 1/1000 of a cell: the corner coordinates in GRIB2 are rounded to
    1e-6 degrees, so computed cell edges can be ~3e-7 degrees off the nominal ones and a box
    drawn exactly on a cell edge must not pull in the neighbouring row or column.
    """
    west, south, east, north = bbox
    if west > east or south > north:
        raise ConversionError("invalid extent (min > max)")
    lats, lons = grid.lats(), grid.lons()
    dy = (lats[0] - lats[-1]) / (len(lats) - 1)
    dx = (lons[-1] - lons[0]) / (len(lons) - 1)
    if tol is None:
        tol = 1e-3 * min(dx, dy)
    # a zero-width/height box (point feature, vertical/horizontal line) takes the cell(s)
    # touching it, including both neighbours when it lies exactly on a cell edge
    tx = -tol if east - west <= 2 * tol else tol
    ty = -tol if north - south <= 2 * tol else tol
    cols = np.flatnonzero((lons + dx / 2 > west + tx) & (lons - dx / 2 < east - tx))
    rows = np.flatnonzero((lats + dy / 2 > south + ty) & (lats - dy / 2 < north - ty))
    if cols.size == 0 or rows.size == 0:
        raise ConversionError(
            f"extent {west:.4f},{south:.4f} - {east:.4f},{north:.4f} does not overlap the grid "
            f"({lons[0] - dx / 2:.4f},{lats[-1] - dy / 2:.4f} - "
            f"{lons[-1] + dx / 2:.4f},{lats[0] + dy / 2:.4f})")
    clipped = (west < lons[0] - dx / 2 - tol or east > lons[-1] + dx / 2 + tol
               or south < lats[-1] - dy / 2 - tol or north > lats[0] + dy / 2 + tol)
    j0, j1 = _at_least_two(int(rows[0]), int(rows[-1]) + 1, len(lats))
    i0, i1 = _at_least_two(int(cols[0]), int(cols[-1]) + 1, len(lons))
    return j0, j1, i0, i1, clipped


def _at_least_two(a: int, b: int, n: int):
    """GDAL cannot georeference a netCDF raster that is 1 cell wide or high: add a neighbour."""
    if b - a >= 2 or n < 2:
        return a, b
    return (a, b + 1) if b < n else (a - 1, b)


def _attr(obj, name, value):
    if isinstance(value, str):
        a = obj.CreateAttribute(name, [], gdal.ExtendedDataType.CreateString())
    else:
        dt = gdal.GDT_Float64 if isinstance(value, float) else gdal.GDT_Int32
        a = obj.CreateAttribute(name, [], gdal.ExtendedDataType.Create(dt))
    a.Write(value)


def _crs_wkt(name, a, rf, epsg):
    return (f'GEOGCRS["{name}",DATUM["{name}",ELLIPSOID["{name} ellipsoid",{a},{rf}]],'
            'CS[ellipsoidal,2],AXIS["latitude",north],AXIS["longitude",east],'
            f'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",{epsg}]]')


def _write(out_path, grid, variables, times, zlevel, reporter, n_files, win=None):
    gdal.UseExceptions()
    drv = gdal.GetDriverByName("netCDF")
    if drv is None:
        raise ConversionError("GDAL netCDF driver is not available")
    ds = drv.CreateMultiDimensional(out_path, [], ["FORMAT=NC4"])
    rg = ds.GetRootGroup()
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    _attr(rg, "Conventions", "CF-1.8")
    _attr(rg, "title", "JMA GRIB2 (run-length packed) converted to NetCDF")
    _attr(rg, "institution", "Japan Meteorological Agency (original data)")
    _attr(rg, "source", f"{n_files} GRIB2 file(s), data representation template 5.200")
    _attr(rg, "history", f"{now} created by {TOOL_NAME}")
    _attr(rg, "time_coverage_start", times[0].strftime("%Y-%m-%dT%H:%M:%SZ"))
    _attr(rg, "time_coverage_end", times[-1].strftime("%Y-%m-%dT%H:%M:%SZ"))

    f64 = gdal.ExtendedDataType.Create(gdal.GDT_Float64)
    f32 = gdal.ExtendedDataType.Create(gdal.GDT_Float32)
    j0, j1, i0, i1 = win if win else (0, grid.nj, 0, grid.ni)
    nt, nj, ni = len(times), j1 - j0, i1 - i0
    lats, lons = grid.lats()[j0:j1], grid.lons()[i0:i1]
    hy = (grid.lats()[0] - grid.lats()[-1]) / (grid.nj - 1) / 2
    hx = (grid.lons()[-1] - grid.lons()[0]) / (grid.ni - 1) / 2
    _attr(rg, "geospatial_lat_min", float(lats[-1] - hy))
    _attr(rg, "geospatial_lat_max", float(lats[0] + hy))
    _attr(rg, "geospatial_lon_min", float(lons[0] - hx))
    _attr(rg, "geospatial_lon_max", float(lons[-1] + hx))
    if win:
        _attr(rg, "comment", f"spatial subset of the {grid.ni} x {grid.nj} source grid "
                             f"(columns {i0}-{i1 - 1}, rows {j0}-{j1 - 1} from the north-west)")
    d_t = rg.CreateDimension("time", "TEMPORAL", None, nt)
    d_y = rg.CreateDimension("lat", "HORIZONTAL_Y", "NORTH", nj)
    d_x = rg.CreateDimension("lon", "HORIZONTAL_X", "EAST", ni)

    v_t = rg.CreateMDArray("time", [d_t], f64)
    v_t.Write(np.array([(t - EPOCH).total_seconds() / 60.0 for t in times]))
    _attr(v_t, "standard_name", "time")
    _attr(v_t, "long_name", "valid time (UTC)")
    _attr(v_t, "units", "minutes since 1970-01-01 00:00:00")
    _attr(v_t, "calendar", "standard")
    _attr(v_t, "axis", "T")

    intervals = {v.key[4] for v in variables.values()}
    if len(intervals) == 1 and None not in intervals:
        iv = intervals.pop()
        d_nv = rg.CreateDimension("nv", None, None, 2)
        v_b = rg.CreateMDArray("time_bnds", [d_t, d_nv], f64)
        tv = np.array([(t - EPOCH).total_seconds() / 60.0 for t in times])
        v_b.Write(np.stack([tv - iv, tv], axis=1))
        _attr(v_t, "bounds", "time_bnds")

    v_y = rg.CreateMDArray("lat", [d_y], f64)
    v_y.Write(lats)
    _attr(v_y, "standard_name", "latitude")
    _attr(v_y, "units", "degrees_north")
    _attr(v_y, "axis", "Y")
    v_x = rg.CreateMDArray("lon", [d_x], f64)
    v_x.Write(lons)
    _attr(v_x, "standard_name", "longitude")
    _attr(v_x, "units", "degrees_east")
    _attr(v_x, "axis", "X")

    a, rf, cname, epsg = _EARTH.get(grid.earth_shape, _EARTH[4])
    if grid.earth_shape not in _EARTH:
        reporter.warn(f"earth shape code {grid.earth_shape} unknown; GRS80/JGD2000 assumed")
    v_c = rg.CreateMDArray("crs", [], gdal.ExtendedDataType.Create(gdal.GDT_Int32))
    _attr(v_c, "grid_mapping_name", "latitude_longitude")
    _attr(v_c, "semi_major_axis", a)
    _attr(v_c, "inverse_flattening", rf)
    _attr(v_c, "crs_wkt", _crs_wkt(cname, a, rf, epsg))

    opts = [f"BLOCKSIZE=1,{nj},{ni}"]
    if zlevel > 0:
        opts += ["COMPRESS=DEFLATE", f"ZLEVEL={int(zlevel)}"]
    arrays = []
    for var in variables.values():
        arr = rg.CreateMDArray(var.name, [d_t, d_y, d_x], f32, opts)
        arr.SetNoDataValueDouble(FILL_VALUE)
        for k, v in var.attrs.items():
            if k == "units":
                arr.SetUnit(v)
            else:
                _attr(arr, k, v)
        _attr(arr, "grid_mapping", "crs")
        _attr(arr, "coordinates", "time lat lon")
        arrays.append((var, arr))

    total = nt * len(arrays)
    done = 0
    empty = np.full((nj, ni), FILL_VALUE, dtype=np.float32)
    for var, arr in arrays:
        for k, t in enumerate(times):
            if reporter.canceled():
                ds = None
                return False
            m = var.messages.get(t)
            if m is None:
                block = empty
            else:
                try:
                    d = decode(m)[j0:j1, i0:i1]
                    block = np.where(np.isnan(d), FILL_VALUE, d).astype(np.float32)
                except (Grib2Error, OSError, ValueError) as e:
                    reporter.warn(f"{os.path.basename(m.path)}: decode failed ({e}); "
                                  "time step written as missing")
                    block = empty
            arr.Write(block, array_start_idx=[k, 0, 0], count=[1, nj, ni])
            done += 1
            reporter.progress(10.0 + 90.0 * done / total)
    ds = None
    return True


def convert(paths: Sequence[str], out_path: str, reporter: Optional[Reporter] = None,
            start: Optional[datetime] = None, end: Optional[datetime] = None,
            zlevel: int = 4, bbox: Optional[Sequence[float]] = None) -> dict:
    """Convert GRIB2 files to one NetCDF file.

    Times are naive datetimes in UTC. bbox = (west, south, east, north) in degrees of the
    grid's geographic CRS limits the output to the grid cells intersecting that box.
    """
    reporter = reporter or Reporter()
    if not paths:
        raise ConversionError("no input files")
    reporter.info(f"scanning {len(paths)} file(s)")
    grid, variables, n_skipped = collect(paths, reporter, start, end)
    if reporter.canceled():
        return dict(canceled=True)
    if grid is None or not any(v.messages for v in variables.values()):
        raise ConversionError("no supported message found in the selected period")
    times = sorted({t for v in variables.values() for t in v.messages})
    for v in variables.values():
        reporter.info(f"variable '{v.name}': {len(v.messages)} time step(s), "
                      f"{v.attrs['long_name']} [{v.attrs['units']}]")
    reporter.info(f"time: {times[0]:%Y-%m-%d %H:%M} .. {times[-1]:%Y-%m-%d %H:%M} UTC, "
                  f"{len(times)} step(s)")
    win = None
    if bbox is not None:
        j0, j1, i0, i1, clipped = subset_window(grid, bbox)
        win = (j0, j1, i0, i1)
        if clipped:
            reporter.warn("the extent reaches beyond the grid; "
                          "only the overlapping part is written")
        lats, lons = grid.lats(), grid.lons()
        reporter.info(f"spatial subset: {i1 - i0} x {j1 - j0} cells "
                      f"(lon {lons[i0]:.5f}..{lons[i1 - 1]:.5f}, "
                      f"lat {lats[j1 - 1]:.5f}..{lats[j0]:.5f}; cell centres), "
                      f"{100.0 * (i1 - i0) * (j1 - j0) / (grid.ni * grid.nj):.2f}% of the grid")

    out_dir = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(out_dir, exist_ok=True)
    if os.path.exists(out_path):
        os.remove(out_path)
    n_files = len(paths) - n_skipped
    target, tmp_dir = out_path, None
    try:
        try:
            ok = _write(target, grid, variables, times, zlevel, reporter, n_files, win)
        except RuntimeError as e:
            # netCDF-C on Windows may fail on non-ASCII paths: retry via a temp file
            reporter.warn(f"direct write failed ({e}); retrying through a temporary file")
            if os.path.exists(out_path):
                os.remove(out_path)
            tmp_dir = tempfile.mkdtemp(prefix="grib2nc_")
            target = os.path.join(tmp_dir, "out.nc")
            ok = _write(target, grid, variables, times, zlevel, reporter, n_files, win)
        if not ok:
            if os.path.exists(target):
                os.remove(target)
            return dict(canceled=True)
        if target != out_path:
            shutil.move(target, out_path)
    finally:
        if tmp_dir:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    size_mb = os.path.getsize(out_path) / 1e6
    reporter.info(f"written {out_path} ({size_mb:.1f} MB)")
    return dict(canceled=False, output=out_path, variables=[v.name for v in variables.values()],
                times=len(times), skipped=n_skipped, window=win)
