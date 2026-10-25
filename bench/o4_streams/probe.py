"""Bench (PLAN D6): transmit-antenna statistics (digital_vtx.streams) on the DJI O4 video bursts of the Mini 4 Pro
testdata capture (D_20251119-172214; Mini 4 Pro = 2T4R per DJI)."""
import glob, json, sys
import numpy as np
sys.path.insert(0, "python")
from spectrum_intel.core import common as C
from spectrum_intel.receivers.dji import detector as DJ
from spectrum_intel.receivers.digital_vtx import streams as ST, detector as DV

meta = sorted(glob.glob("testdata/DJI_Mini_Pro_4_20251119-172214/*.sigmf-meta"))[0]
x, fs, fc, _ = C.load_capture(meta)
x = x[:1 << 21]
em = [e for e in DJ.detect(x, fs, fc) if e.label.startswith("video_o4")]
print(len(em), "O4 bursts", em and em[0].features.get("carrier_mhz"), em and em[0].features.get("bw_mhz"))
f = em[0].features
X = np.fft.fft(x.astype(np.complex128))
tu, cp = f["tu_us"] * 1e-6, f["cp_us"] * 1e-6
ofs = 1024 / tu                                                    # 30.72 MSps: Tu = 1024 samples
y, ofs = DV._slice(X, fs, f["carrier_mhz"] * 1e6 - fc, f["bw_mhz"] * 1e6 / 2 + 0.3e6, ofs)
runs = [(e.start_us * 1e-6, e.end_us * 1e-6, e.features.get("n_symbols", 0)) for e in em]
out = ST.tx_streams(y, ofs, tu, tu + cp, f["bw_mhz"] * 1e6, runs, 20.0, len(x) / fs)
print(json.dumps(out, indent=1, default=float)[:1500])
