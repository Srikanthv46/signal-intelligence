"""Bench (PLAN D2): frame starts of the Walksnail-like link in LIVE_OFDM_VTX_5805MHz (digital_vtx emissions)."""
import sys
sys.path.insert(0, "python")
from spectrum_intel.core import common as C
from spectrum_intel.receivers.digital_vtx import detector as D

x, fs, fc, _ = C.load_capture("testdata/live/LIVE_OFDM_VTX_5805MHz/LIVE_OFDM_VTX_5805MHz.sigmf-meta")
x = x[:1 << 21]
em = [e for e in D.detect(x, fs, fc) if e.label.startswith("ofdm_link")]
for e in em:
    print(round(e.start_us, 1), round(e.end_us, 1), e.center_mhz, {k: e.features.get(k) for k in ("tu_us", "cp_us", "carrier_mhz", "frame_period_ms", "bw_mhz")})
