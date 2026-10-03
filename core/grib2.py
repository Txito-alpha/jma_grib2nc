# -*- coding: utf-8 -*-
"""Minimal GRIB2 reader for JMA run-length packed products (data template 5.200).

QGIS independent: depends on numpy only.

GDAL (all versions) cannot decode template 5.200, and ecCodes cannot parse JMA's
local product definition template 4.50008, so the sections needed here are
parsed directly.
"""
import struct
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple

import numpy as np

SUPPORTED_GRID_TEMPLATES = (0,)
SUPPORTED_DATA_TEMPLATES = (200,)

# Code table 4.4 (indicator of unit of time range) -> seconds
_TIME_UNIT_SECONDS = {
    0: 60, 1: 3600, 2: 86400, 10: 3 * 3600, 11: 6 * 3600, 12: 12 * 3600, 13: 1,
}
# Product definition templates whose octets 35-41 hold the end of the time interval
_PDT_WITH_INTERVAL_END = (8, 50008)


class Grib2Error(Exception):
    """Raised when a message cannot be read by this reader."""


@dataclass(frozen=True)
class Grid:
    ni: int
    nj: int
    lat1: float
    lon1: float
    lat2: float
    lon2: float
    dx: float
    dy: float
    scan: int
    earth_shape: int

    def lats(self) -> np.ndarray:
        """Row latitudes in output order (north to south).

        The increments stored in the file are rounded to 1e-6 degrees
        (1/120 deg -> 0.008333), so coordinates are derived from the corner
        points instead of from the increment (avoids ~120 m drift at 3360 rows).
        """
        north, south = max(self.lat1, self.lat2), min(self.lat1, self.lat2)
        return np.linspace(north, south, self.nj)

    def lons(self) -> np.ndarray:
        return np.linspace(min(self.lon1, self.lon2), max(self.lon1, self.lon2), self.ni)


@dataclass
class MessageInfo:
    path: str
    offset: int
    length: int
    discipline: int
    reference_time: datetime
    valid_time: datetime
    grid: Grid
    pdt: int
    drt: int
    category: int
    number: int
    stat_process: Optional[int] = None
    interval_minutes: Optional[int] = None
    extra: Dict[str, object] = field(default_factory=dict)

    @property
    def param_key(self) -> Tuple[int, int, int, Optional[int], Optional[int]]:
        return (self.discipline, self.category, self.number,
                self.stat_process, self.interval_minutes)


def _s(b: bytes) -> int:
    """GRIB2 signed integer (sign bit + magnitude, not two's complement)."""
    v = int.from_bytes(b, "big")
    top = 1 << (len(b) * 8 - 1)
    return -(v - top) if v & top else v


def _u(b: bytes) -> int:
    return int.from_bytes(b, "big")


def _dt(b: bytes) -> datetime:
    y, mo, d, h, mi, s = struct.unpack(">HBBBBB", b[:7])
    return datetime(y, mo, d, h, mi, s)


def _split_sections(buf: bytes, start: int) -> Tuple[Dict[int, bytes], int]:
    if buf[start:start + 4] != b"GRIB":
        raise Grib2Error("GRIB indicator not found")
    edition = buf[start + 7]
    if edition != 2:
        raise Grib2Error(f"GRIB edition {edition} is not supported (GRIB2 only)")
    total = _u(buf[start + 8:start + 16])
    end = start + total
    if end > len(buf):
        raise Grib2Error("truncated message")
    if buf[end - 4:end] != b"7777":
        raise Grib2Error("end section '7777' not found")
    secs = {0: buf[start:start + 16]}
    p = start + 16
    while p < end - 4:
        length = _u(buf[p:p + 4])
        num = buf[p + 4]
        if length < 5 or p + length > end:
            raise Grib2Error(f"broken section {num} at byte {p}")
        secs[num] = buf[p:p + length]
        p += length
    for need in (1, 3, 4, 5, 7):
        if need not in secs:
            raise Grib2Error(f"section {need} missing")
    return secs, end


def parse_grid(s3: bytes) -> Grid:
    tmpl = _u(s3[12:14])
    if tmpl not in SUPPORTED_GRID_TEMPLATES:
        raise Grib2Error(f"grid definition template 3.{tmpl} is not supported")
    return Grid(
        ni=_u(s3[30:34]), nj=_u(s3[34:38]),
        lat1=_s(s3[46:50]) / 1e6, lon1=_s(s3[50:54]) / 1e6,
        lat2=_s(s3[55:59]) / 1e6, lon2=_s(s3[59:63]) / 1e6,
        dx=_u(s3[63:67]) / 1e6, dy=_u(s3[67:71]) / 1e6,
        scan=s3[71], earth_shape=s3[14],
    )


def _parse_product(s4: bytes, ref: datetime) -> dict:
    pdt = _u(s4[7:9])
    out = dict(pdt=pdt, category=s4[9], number=s4[10])
    unit = _TIME_UNIT_SECONDS.get(s4[17])
    ft = _s(s4[18:22])
    if pdt in _PDT_WITH_INTERVAL_END and len(s4) >= 53:
        # end of overall time interval = valid time of an accumulated product
        out["valid_time"] = _dt(s4[34:41])
        out["stat_process"] = s4[46]
        iv_unit = _TIME_UNIT_SECONDS.get(s4[48])
        if iv_unit:
            out["interval_minutes"] = _u(s4[49:53]) * iv_unit // 60
    elif unit is not None:
        out["valid_time"] = ref + timedelta(seconds=ft * unit)
    else:
        out["valid_time"] = ref
        out["warning"] = f"unknown time unit {s4[17]}; reference time used as valid time"
    return out


def _scan_one(buf: bytes, start: int, path: str) -> Tuple[MessageInfo, int]:
    secs, end = _split_sections(buf, start)
    ref = _dt(secs[1][12:19])
    prod = _parse_product(secs[4], ref)
    info = MessageInfo(
        path=path, offset=start, length=end - start,
        discipline=secs[0][6], reference_time=ref, valid_time=prod["valid_time"],
        grid=parse_grid(secs[3]), pdt=prod["pdt"], drt=_u(secs[5][9:11]),
        category=prod["category"], number=prod["number"],
        stat_process=prod.get("stat_process"),
        interval_minutes=prod.get("interval_minutes"),
    )
    if "warning" in prod:
        info.extra["warning"] = prod["warning"]
    return info, end


def scan_file(path: str) -> List[MessageInfo]:
    """Read headers of all GRIB2 messages in a file (data are not decoded)."""
    with open(path, "rb") as f:
        buf = f.read()
    msgs, p = [], 0
    while True:
        start = buf.find(b"GRIB", p)
        if start < 0:
            break
        info, p = _scan_one(buf, start, path)
        msgs.append(info)
    if not msgs:
        raise Grib2Error("no GRIB message found")
    return msgs


def decode_run_length(s5: bytes, s7: bytes, npts: int) -> np.ndarray:
    """Decode data template 5.200 (run length packing with level values).

    Returns float64 values; level 0 (missing) becomes NaN.
    """
    tmpl = _u(s5[9:11])
    if tmpl != 200:
        raise Grib2Error(f"data representation template 5.{tmpl} is not supported")
    nbit = s5[11]
    vmax = _u(s5[12:14])
    mlev = _u(s5[14:16])
    dscale = _s(s5[16:17])
    if len(s5) < 17 + 2 * mlev:
        raise Grib2Error("section 5 shorter than its level table")
    levels = np.array([_s(s5[17 + 2 * i:19 + 2 * i]) for i in range(mlev)],
                      dtype=np.float64) / (10.0 ** dscale)

    raw = np.frombuffer(s7[5:], dtype=np.uint8)
    if nbit == 8:
        vals = raw.astype(np.int64)
    elif nbit == 16:
        vals = raw[: len(raw) // 2 * 2].view(">u2").astype(np.int64)
    else:
        bits = np.unpackbits(raw)
        n = len(bits) // nbit
        weights = 1 << np.arange(nbit - 1, -1, -1, dtype=np.int64)
        vals = bits[: n * nbit].reshape(n, nbit).astype(np.int64).dot(weights)
    lngu = (1 << nbit) - 1 - vmax
    if lngu <= 0:
        raise Grib2Error("invalid run length parameters")

    is_lv = vals <= vmax
    lv_idx = np.flatnonzero(is_lv)
    if lv_idx.size == 0 or lv_idx[0] != 0:
        raise Grib2Error("run length stream does not start with a level value")
    grp = np.cumsum(is_lv) - 1                     # level-value group of each item
    pos = np.arange(len(vals)) - lv_idx[grp] - 1   # digit position within a run
    digit = np.where(is_lv, 0, vals - (vmax + 1))
    contrib = digit * np.power(lngu, np.maximum(pos, 0), dtype=np.int64)
    runs = np.bincount(grp, weights=contrib, minlength=lv_idx.size).astype(np.int64) + 1

    # trailing padding bits can form spurious items: trim to npts
    counts = np.cumsum(runs)
    if counts[-1] < npts:
        raise Grib2Error(f"decoded {counts[-1]} points, expected {npts}")
    out_lv = np.repeat(vals[lv_idx], runs)[:npts]
    if np.any(out_lv > mlev):
        raise Grib2Error("level value exceeds level table")
    lut = np.concatenate([[np.nan], levels])
    return lut[out_lv]


def decode(info: MessageInfo) -> np.ndarray:
    """Decode one message to a (nj, ni) float64 array, north row first."""
    with open(info.path, "rb") as f:
        f.seek(info.offset)
        buf = f.read(info.length)
    secs, _ = _split_sections(buf, 0)
    if 6 in secs and secs[6][5] != 255:
        raise Grib2Error("bit-map section is not supported")
    g = info.grid
    data = decode_run_length(secs[5], secs[7], g.ni * g.nj)
    if g.scan & 0x20:
        raise Grib2Error("scanning mode with adjacent points in j direction is not supported")
    arr = data.reshape(g.nj, g.ni)
    if g.scan & 0x40:   # rows run south -> north
        arr = arr[::-1]
    if g.scan & 0x80:   # columns run east -> west
        arr = arr[:, ::-1]
    return arr
