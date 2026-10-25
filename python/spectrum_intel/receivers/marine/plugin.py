"""Pipeline plugin: AIS vessel / station reports (wraps .detector; research, signal gaps S15)."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver

MERGE = ("flag", "mid", "kind", "nav_status", "sog_kn", "lat", "lon", "cog_deg", "heading_deg", "name", "callsign",
         "ship_type", "destination", "imo", "length_m", "beam_m", "draught_m", "aid_type", "utc")


class MarineDecoder(EmissionReceiver):
    name = "marine"
    protocol = "AIS"
    key = "marine"
    bands = {"marine_vhf"}
    method = "receivers.marine: AIS GMSK / HDLC / CRC-16 X.25 (research)"
    owns_band = True
    steps = ["AIS 1 / 2 (161.975 / 162.025 MHz, receiver ppm applied): +-10 kHz at 48 kSps, FM discriminator, "
             "offset removed (32-bit moving average), integrated over one bit",
             "per sample phase: NRZ levels -> NRZI decoding -> training + flag 0x7E, de-stuffing to the closing "
             "flag, octets LSB first, CRC-16/X.25 residue 0xF0B8",
             "message fields (ITU-R M.1371-5): types 1/2/3 position, 4 base station, 5 static / voyage, 18 class B, "
             "21 aid to navigation, 24 class B static; MID 419 = India"]

    def run_detector(self, ctx):
        from .detector import detect
        cal = ctx.cache.get("ppm") or ctx.cache.get("sweep_ppm")
        return detect(ctx.ch, cal["ppm"] if cal else 0.0)

    def group_key(self, e):
        return "ais"

    def summarise(self, em, sig):
        vessels: dict = {}
        msgs = []
        for e in sorted(em, key=lambda e: e.start_us):
            d = e.decoded
            msgs.append({"t_s": round(e.start_us * 1e-6, 5), "channel": e.features["channel"],
                         **{k: v for k, v in d.items() if k != "crc_ok"}})
            v = vessels.setdefault(d["mmsi"], {"mmsi": d["mmsi"], "messages": 0, "msg_types": []})
            v["messages"] += 1
            if d["msg_type"] not in v["msg_types"]:
                v["msg_types"].append(d["msg_type"])
            for k in MERGE:
                if d.get(k) is not None:
                    v[k] = d[k]
        vs = list(vessels.values())
        names = [v.get("name") or v["mmsi"] for v in vs]
        fields = {"n_messages": len(msgs), "vessels": vs, "messages": msgs[:100]}
        return result("success", self.method, self.steps, fields, "every message verified by its CRC-16/X.25",
                      label=f"AIS ({len(vs)} station{'s' if len(vs) != 1 else ''})", confidence=0.97,
                      modulation="GMSK 9600 bit/s (AIS)",
                      evidence=f"{len(msgs)} CRC-valid AIS message(s) from {names}",
                      span=(min(e.f_lo_mhz for e in em) * 1e6, max(e.f_hi_mhz for e in em) * 1e6))


PLUGINS = [MarineDecoder()]
