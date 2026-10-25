"""Pipeline plugin: FM broadcast stations and RDS (wraps .detector; research, signal gaps S20)."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean


class BroadcastFmDecoder(EmissionReceiver):
    name = "broadcast_fm"
    protocol = "FM broadcast + RDS"
    key = "broadcast_fm"
    bands = {"fm_bcast"}
    method = "receivers.broadcast_fm: FM discriminator, stereo pilot, RDS block code (research)"
    owns_band = True                 # the station's sidebands would otherwise become 'unknown' boxes
    steps = ["100 kHz raster 87.6-107.9 MHz: channel power over +-100 kHz >= 10 dB over the band floor (local "
             "maxima on the raster)",
             "+-120 kHz at 240 kSps, FM discriminator -> MPX: peak deviation 15-110 kHz; stereo when the 19 kHz "
             "pilot line is >= 10 dB over its surroundings",
             "RDS (blocks >= 45 ms): MPX shifted by 3 x the pilot, +-2.4 kHz, BPSK phase from the squared signal, "
             "biphase matched filter at 1187.5 bit/s, differential decoding, block sync on the offset words "
             "(error-free (26,16) blocks only) -> PI, PTY, TP, programme name, RadioText, clock time"]

    def run_detector(self, ctx):
        from .detector import detect
        return detect(ctx.ch)

    def group_key(self, e):
        return e.features["channel_mhz"]

    def summarise(self, em, sig):
        e = max(em, key=lambda e: e.confidence)
        F = e.features
        fields = _clean(dict(F))
        rds = F.get("rds") or {}
        parts = [f"{F['channel_mhz']:.1f} MHz"] + (["stereo"] if F["stereo"] else [])
        if rds.get("ps"):
            parts.append(f"RDS '{rds['ps'].replace('_', ' ').strip()}'")
        elif rds.get("pi"):
            parts.append(f"RDS PI {rds['pi']}")
        ok = bool(rds.get("pi"))
        ev = (f"FM carrier {F['level_db']} dB over the band floor, deviation {F['deviation_khz']} kHz"
              + (f", 19 kHz pilot +{F['pilot_db']} dB" if F["stereo"] else ", no stereo pilot")
              + (f"; RDS: {rds['blocks_ok']} error-free block(s), PI {rds['pi']}"
                 + (f", PTY {rds['pty_name']}" if rds.get("pty_name") else "") if ok else ""))
        return result("success" if ok else "partial", self.method, self.steps, fields,
                      "RDS blocks verified by their (26,16) checkwords" if ok else
                      "FM broadcast measured; no RDS in this block (RDS needs >= 45 ms, many stations send none)",
                      label=f"FM broadcast ({', '.join(parts)})", confidence=e.confidence,
                      modulation="wideband FM (+-75 kHz), stereo MPX" if F["stereo"] else "wideband FM (+-75 kHz)",
                      evidence=ev, span=(min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6))


PLUGINS = [BroadcastFmDecoder()]
