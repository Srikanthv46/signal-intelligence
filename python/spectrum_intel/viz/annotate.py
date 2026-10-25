"""Annotated spectrogram drawn from a report JSON (so sweep-level passes can
patch a report and redraw without re-analysing): absolute time / frequency
axes, one box per burst, label per signal, colour by decode status."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from ..detect.stft import stft  # noqa: E402

COLORS = {"success": "#00e676", "partial": "#ffd600", "not_possible": "#ff5252", "artifact": "#b0bec5"}
ARTIFACT_FLAGS = {"lo_leakage", "iq_image", "receiver_spur", "lo_harmonic_response"}


def _color(s):
    if ARTIFACT_FLAGS & set(s["flags"]):
        return COLORS["artifact"]
    return COLORS.get(s["decoded"]["status"], "#ffffff")


def _short(s):
    f = s["decoded"]["fields"]
    for k in ("ssids", "pci", "icao", "bsic", "pcis"):
        v = f.get(k)
        if v not in (None, [], ""):
            v = ",".join(map(str, v[:3])) if isinstance(v, list) else v
            return f" {k.rstrip('s').upper() if k != 'ssids' else 'SSID'} {v}"
    return ""


def draw(rep: dict, iq: np.ndarray, path: Path, max_boxes_per_signal: int = 60):
    fs, fc = rep["sample_rate_hz"], rep["center_freq_hz"]
    C = stft(iq, fs, fc, 1024)
    db = 10 * np.log10(C.P.T + 1e-30)
    fl = float(np.percentile(db, 20))
    fig, ax = plt.subplots(figsize=(16, 9.5))
    ax.imshow(db, origin="lower", aspect="auto", cmap="viridis", vmin=fl - 2, vmax=fl + 40,
              extent=[0, (C.t[-1] + C.frame_len) * 1e3, C.f[0] / 1e6, C.f[-1] / 1e6], interpolation="nearest")
    for s in rep["signals"]:
        col = _color(s)
        bs = s["features"].get("bursts") or [{"t_start_s": s["t_start_s"], "t_end_s": s["t_end_s"],
                                              "f_low_hz": s["f_low_hz"], "f_high_hz": s["f_high_hz"]}]
        ls = "-"
        if len(bs) > max_boxes_per_signal:
            bs, ls = [{"t_start_s": s["t_start_s"], "t_end_s": s["t_end_s"],
                       "f_low_hz": s["f_low_hz"], "f_high_hz": s["f_high_hz"]}], "--"
        for b in bs:
            ax.add_patch(Rectangle((b["t_start_s"] * 1e3, b["f_low_hz"] / 1e6),
                                   (b["t_end_s"] - b["t_start_s"]) * 1e3, (b["f_high_hz"] - b["f_low_hz"]) / 1e6,
                                   fill=False, ec=col, lw=1.0, ls=ls))
        lab = f"#{s['id']} {s['label']} ({s['confidence']:.2f}){_short(s)}"
        ax.text(s["t_start_s"] * 1e3 + 0.1, s["f_high_hz"] / 1e6, lab, color=col, fontsize=7, va="bottom",
                bbox=dict(facecolor="black", alpha=0.55, pad=1, lw=0))
    ax.set_xlabel("time since capture start (ms)")
    ax.set_ylabel("frequency (MHz)")
    lc = rep.get("frequency_label_check", {})
    ttl = (f"{rep['file']}   fc {fc / 1e6:.2f} MHz   fs {fs / 1e6:.2f} MSps   {rep['capture_time_utc'][:23]}   "
           f"{len(rep['signals'])} signals")
    if str(lc.get("status", "")).startswith("corrected"):
        ttl += (f"\nNOTE: metadata fc {lc['fc_metadata_hz'] / 1e6:.2f} MHz is wrong ({lc['offset_mhz']:+.1f} MHz "
                f"from {lc.get('n_adverts')} BLE adverts); IQ analysed at {fc / 1e6:.2f} MHz")
    q = rep.get("capture_quality", {})
    if q.get("clipping"):
        ttl += f"   ADC CLIPPING ({q['clipped_fraction']:.1e} of samples)"
    ax.set_title(ttl, fontsize=10)
    hs = [Rectangle((0, 0), 1, 1, fill=False, ec=c) for c in COLORS.values()]
    ax.legend(hs, ["decoded", "partially decoded", "not decoded", "receiver artifact"], loc="upper right",
              fontsize=8, framealpha=0.6)
    fig.tight_layout()
    # 8-bit palette PNG: several times smaller than RGBA for a noise-textured spectrogram
    import io
    from PIL import Image
    buf = io.BytesIO()
    fig.savefig(buf, dpi=90, format="png")
    plt.close(fig)
    buf.seek(0)
    Image.open(buf).convert("RGB").quantize(colors=128, method=Image.Quantize.MEDIANCUT).save(path, optimize=True)
    return path
