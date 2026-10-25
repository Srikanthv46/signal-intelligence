"""Pipeline plugin: DVB-T-framed COFDM links (drone / UAV video downlinks; research capability) (wraps .detector)."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean


class CofdmDecoder(EmissionReceiver):
    name = "cofdm"
    protocol = "DVB-T COFDM video link"
    key = "cofdm"
    method = "receivers.cofdm (research receiver: ETSI EN 300 744, synthetic validation only)"
    steps = ["candidates: Welch spectrum runs 0.951 x B wide (B = 8 / 7 / 6 / 5 / 4 / 2.5 / 2 MHz) above floor + 2 dB",
             "numerology: slice at 8 B / 7 MHz, CP metric folded on N (1 + g) for 2K / 8K x GI 1/4..1/32",
             "carrier: fractional offset from the CP phase, integer offset from the 4/3-boosted continual pilots",
             "TPS: DBPSK across symbols on the 17 / 68 TPS carriers, sync word, BCH(67,53) -> mode, GI, "
             "constellation, hierarchy, code rates, frame number, cell-ID byte"]
    MIN_RATE_HZ = 2.5e6

    def run_detector(self, ctx):
        from .detector import detect
        if min(ctx.cap.fs, ctx.cap.bandwidth_hz or ctx.cap.fs) < self.MIN_RATE_HZ:
            return []
        return detect(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc)

    def in_range(self, lo, hi):
        return True                      # COFDM downlinks: 200 MHz - 6 GHz, any band

    def group_key(self, e):
        return e.center_mhz

    def summarise(self, em, sig):
        e = em[0]
        F = e.features
        if F.get("p1"):
            # DVB-T2 P1 preamble (research, PLAN D5): S1 / S2 decoded by pattern correlation (no CRC)
            p1 = F["p1"]
            keys = ("standard", "channel_bw_mhz", "occupied_mhz", "elementary_rate_mhz", "carrier_mhz", "p1_count",
                    "p1_metric", "p1_start_us", "p1")
            fields = _clean({k: F[k] for k in keys if k in F})
            B = f"{F['channel_bw_mhz']:g} MHz"
            like = "" if e.label == "dvbt2" else "-like"
            label = f"COFDM video link{like} (DVB-T2 {p1['preamble']}, {B}, {p1['fft_gi']})"
            ev = (f"P1 C-A-B correlation {F['p1_metric']} ({F['p1_count']} in the block); S1 {p1['s1']} / S2 {p1['s2']} "
                  f"(pattern margins {p1['s1_margin']} / {p1['s2_margin']})")
            return result("success" if not like else "partial", self.method, self.steps, fields,
                          "P1 signalling decoded (S1 / S2 by pattern correlation); L1 / payload not decoded",
                          label=label, confidence=e.confidence, modulation="COFDM (DVB-T2)", evidence=ev,
                          span=(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6))
        tps = F.get("tps") or {}
        keys = ("standard", "channel_bw_mhz", "occupied_mhz", "mode", "guard_interval", "fft_size", "elementary_rate_mhz",
                "scs_hz", "tu_us", "symbol_us", "carrier_mhz", "integer_offset_carriers", "cp_peak", "cp_contrast",
                "symbols", "continual_pilot_db", "continual_pilot_margin", "tps", "tps_note")
        fields = _clean({k: F[k] for k in keys if k in F})
        B = f"{F['channel_bw_mhz']:g} MHz"
        if tps.get("bch_ok"):
            label = (f"COFDM video link (DVB-T {F['mode']}, {B}, {tps['constellation']}, HP {tps['code_rate_hp']}, "
                     f"GI {tps['guard_interval']})")
            status, reason = "success", "TPS decoded (BCH(67,53) passed); the video payload (MPEG-TS, often AES) is not decoded"
        else:
            label = f"COFDM video link-like (DVB-T {F['mode']} numerology, {B}, GI {F['guard_interval']})"
            status, reason = "partial", F.get("tps_note") or "TPS not BCH-verified"
        ev = (f"CP fold peak {F['cp_peak']} ({F['cp_contrast']}x median) at {F['mode']} GI {F['guard_interval']}, "
              f"{F['symbols']} symbols; continual pilots +{F['continual_pilot_db']} dB at offset {F['integer_offset_carriers']}")
        if tps:
            ev += f"; TPS frame {tps.get('frame_number')}, BCH {'ok' if tps.get('bch_ok') else 'failed'}"
        return result(status, self.method, self.steps, fields, reason, label=label, confidence=e.confidence,
                      modulation=f"COFDM {F['mode']} ({tps.get('constellation', 'QPSK/16/64-QAM')})", evidence=ev,
                      span=(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6))


PLUGINS = [CofdmDecoder()]
