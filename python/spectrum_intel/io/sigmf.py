"""SigMF reader: one Capture per .sigmf-meta/.sigmf-data pair."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

_DTYPES = {
    "cf32_le": np.complex64, "cf64_le": np.complex128,
    "ci16_le": np.int16, "ci8": np.int8, "cu8": np.uint8,
}


@dataclass
class Capture:
    path: Path                      # the .sigmf-meta
    fc: float                       # centre frequency, Hz
    fs: float                       # sample rate, Hz
    t0: datetime | None             # capture start (UTC)
    datatype: str
    meta: dict = field(repr=False)
    _iq: np.ndarray | None = field(default=None, repr=False)
    _bandwidth_hz: float | None = None                # analog bandwidth given by the caller (API)

    @property
    def name(self) -> str:
        return self.path.name[: -len(".sigmf-meta")]

    @property
    def data_path(self) -> Path:
        return self.path.with_name(self.name + ".sigmf-data")

    @property
    def png_path(self) -> Path:
        return self.path.with_name(self.name + ".png")

    @property
    def gain_db(self) -> float | None:
        return self.meta["global"].get("signalvault:gain_db")

    @property
    def bandwidth_hz(self) -> float | None:
        """Analog (front-end) bandwidth, if known: given by the caller, or signal-vault's
        signalvault:rf_bandwidth_hz. None = assume the whole sampled band is usable."""
        if self._bandwidth_hz:
            return self._bandwidth_hz
        bw = self.meta["global"].get("signalvault:rf_bandwidth_hz")
        return float(bw) if bw else None

    @property
    def iq(self) -> np.ndarray:
        """complex64, mean (DC) kept: the LO leakage is reported as a signal."""
        if self._iq is None:
            self._iq = read_iq(self.data_path, self.datatype)
        return self._iq

    @property
    def n(self) -> int:
        return len(self.iq)

    @property
    def duration(self) -> float:
        return self.n / self.fs

    @classmethod
    def from_iq(cls, iq: np.ndarray, fs: float, fc: float, t0: datetime | None = None,
                name: str = "live", gain_db: float | None = None, bandwidth_hz: float | None = None) -> "Capture":
        """An in-memory capture (a live receiver block): same interface as a file."""
        meta = {"global": {"core:datatype": "cf32_le", "core:sample_rate": float(fs),
                           **({"signalvault:gain_db": gain_db} if gain_db is not None else {}),
                           **({"signalvault:rf_bandwidth_hz": float(bandwidth_hz)} if bandwidth_hz else {})},
                "captures": [{"core:frequency": float(fc), "core:sample_start": 0,
                              **({"core:datetime": t0.isoformat()} if t0 else {})}]}
        return cls(path=Path(f"{name}.sigmf-meta"), fc=float(fc), fs=float(fs), t0=t0, datatype="cf32_le",
                   meta=meta, _iq=np.asarray(iq, np.complex64), _bandwidth_hz=bandwidth_hz)


def read_iq(path: Path, datatype: str) -> np.ndarray:
    dt = _DTYPES.get(datatype)
    if dt is None:
        raise ValueError(f"unsupported SigMF datatype {datatype!r}")
    raw = np.fromfile(path, dtype=dt)
    if np.iscomplexobj(raw):
        return raw.astype(np.complex64, copy=False)
    raw = raw.astype(np.float32)
    if datatype == "cu8":
        raw -= 127.5
    return (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)


def load(path: str | Path) -> Capture:
    p = Path(path)
    if p.suffix == ".sigmf-data" or p.suffix == ".png":
        p = p.with_suffix(".sigmf-meta")
    if p.is_dir():
        metas = sorted(p.glob("*.sigmf-meta"))
        if len(metas) != 1:
            raise ValueError(f"{p}: expected one .sigmf-meta, found {len(metas)}")
        p = metas[0]
    meta = json.loads(p.read_text())
    g, cap = meta["global"], meta["captures"][0]
    t0 = cap.get("core:datetime")
    return Capture(
        path=p, fc=float(cap["core:frequency"]), fs=float(g["core:sample_rate"]),
        t0=datetime.fromisoformat(t0.replace("Z", "+00:00")) if t0 else None,
        datatype=g["core:datatype"], meta=meta,
    )


def find_captures(folder: str | Path) -> list[Path]:
    """Every .sigmf-meta under folder, sorted by centre frequency then time."""
    metas = sorted(Path(folder).rglob("*.sigmf-meta"))
    def key(m):
        c = json.loads(m.read_text())["captures"][0]
        return (c["core:frequency"], c.get("core:datetime", ""))
    return sorted(metas, key=key)
