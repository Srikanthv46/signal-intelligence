"""The record every research receiver emits (one burst / one continuous segment),
plus the SigMF loaders they use."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

import numpy as np


# ============================================================== I/O
_DT = {"cf32_le": np.complex64, "cf64_le": np.complex128}

def load_sigmf(meta_path, data_path=None):
    meta = json.load(open(meta_path))
    g, cap = meta["global"], meta["captures"][0]
    data_path = data_path or meta_path.replace(".sigmf-meta", ".sigmf-data")
    dt = g["core:datatype"]
    if dt in _DT:
        x = np.fromfile(data_path, dtype=_DT[dt])
    elif dt in ("ci16_le", "ci8", "ci32_le"):
        base = {"ci16_le": "<i2", "ci8": "i1", "ci32_le": "<i4"}[dt]
        raw = np.fromfile(data_path, dtype=base).astype(np.float32)
        x = (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)
    else:
        raise ValueError(f"unsupported datatype {dt}")
    return x, float(g["core:sample_rate"]), float(cap["core:frequency"]), meta


@dataclass
class Emission:
    """One burst / one continuous transmission segment."""
    family: str                 # dji_ocusync | dji_lightbridge | dji_droneid | herelink | analog_vtx |
                                # digital_vtx | rc_fhss | wifi_drone | unknown_ofdm | ...
    label: str                  # finer label, e.g. "ocusync_video_10MHz", "ocusync_uplink", "droneid"
    start_us: float
    end_us: float
    f_lo_mhz: float
    f_hi_mhz: float
    center_mhz: float
    snr_db: float | None = None
    confidence: float = 0.5     # 0..1
    decoded: dict | None = None # payload if actually demodulated (e.g. DroneID serial/GPS)
    # physical-layer fingerprint used for association (fill what you can measure):
    #   cfo_hz          residual carrier offset vs. the nominal channel centre
    #   cfo_ppm         cfo_hz / carrier (same oscillator => same ppm across hops)
    #   clock_ppm       sample/symbol-clock offset measured from symbol period
    #   symbol_us, cp_us, scs_khz, bw_mhz, burst_us, frame_period_ms, frame_phase_us
    #   zc_root, pilot_id, mac, serial ...
    features: dict = field(default_factory=dict)

    def to_dict(self):
        return asdict(self)

def load_capture(meta_path):
    x, fs, fc, meta = load_sigmf(meta_path)
    return x.astype(np.complex64), fs, fc, meta

def save_emissions(path, emissions, extra=None):
    def clean(o):
        if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [clean(v) for v in o]
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.bool_): return bool(o)
        if isinstance(o, np.ndarray): return clean(o.tolist())
        return o
    json.dump(clean(dict(emissions=[e.to_dict() for e in emissions], **(extra or {}))), open(path, "w"), indent=1)
