"""Pipeline plugin: ADS-B / Mode S replies as capture-level signals."""
from __future__ import annotations

from collections import Counter
import numpy as np

from spectrum_intel.decode.base import Decoder, receiver_signal, result
from . import modeac
from .messages import fields
from .modes import FS, F_ADSB, SPU, demod
from .tracks import tracks


def _receiver_position(cap):
    """(lat, lon[, alt_m]) from SigMF core:geolocation (GeoJSON point), if the capture has one."""
    g = (getattr(cap, "meta", None) or {}).get("global", {}).get("core:geolocation") or {}
    c = g.get("coordinates") if g.get("type") == "Point" else None
    return (c[1], c[0], *c[2:3]) if c and len(c) >= 2 else None


class AdsbDecoder(Decoder):
    name = "adsb"
    protocol = "Mode S / ADS-B 1090 MHz"
    capture_level = True

    def can_decode(self, sig, ctx):
        return sig.f0 - 1e6 <= F_ADSB <= sig.f1 + 1e6

    def capture_signals(self, ctx):
        cap = ctx.cap
        if not (cap.fc - 0.45 * cap.fs < F_ADSB < cap.fc + 0.45 * cap.fs):
            return []
        y, ofs = ctx.ch.slice(F_ADSB, 2e6, FS)
        env = np.abs(y)
        noise = np.median(env)
        thr = 6 * noise * SPU / 2 * 2            # preamble correlation for pulses ~>6x noise amplitude
        msgs = fields(demod(env, thr))
        replies = modeac.decode(env, FS, [(m["t_s"], m["t_s"] + (8 + m["bits"]) * 1e-6) for m in msgs], noise)
        steps = ["channelise 1090 MHz +-2 MHz at 8 MSps, magnitude envelope",
                 "8 us preamble correlation (pulses 0/1/3.5/4.5 us), PPM bit decisions",
                 "CRC-24 (0xFFF409): DF17/18 residual 0; other DFs address/parity",
                 "fields of every DF/typecode (pyModeS + DO-260B bit layouts): identity, category, altitude, "
                 "velocity, intent, integrity, capabilities, Comm-B registers, ACAS and radar reply fields",
                 "per-aircraft: CPR position (even/odd global, then local), BDS50/60 disambiguation, "
                 "wind and air temperature, interrogation environment, plausibility flags"]
        steps_ac = ["SSR Mode A/C replies (research): 0.45 us pulses >= 12 dB over the envelope median, framing pulses "
                    "F1 / F2 20.3 us apart within 3 dB, every pulse between them on the 1.45 us grid (else garbled), "
                    "X empty; SPI at 24.65 us",
                    "each reply read both ways: Mode A identity (octal ABCD) and Mode C Gillham altitude when the "
                    "bits are a valid Gillham code; words repeated in the block are the verified ones"]
        if not msgs:
            ctx.cache["adsb_note"] = "no Mode S frame passed CRC in this capture"
            if not replies:
                return []
            return [self._mode_ac_signal(replies, steps[:1] + steps_ac)]
        icaos = sorted({m["icao"] for m in msgs if m.get("icao")})
        t0 = cap.t0.timestamp() if getattr(cap, "t0", None) else 0.0
        for m in msgs:
            m["t"] = t0 + m["t_s"]
        ref = _receiver_position(cap)
        aircraft = tracks(msgs, ref)
        f = {"n_frames": len(msgs), "dfs": dict(Counter(m["df"] for m in msgs)), "icao": icaos,
             "aircraft": aircraft, "receiver_position": ref, "frames": msgs}
        if replies:
            f["mode_ac"] = modeac.summarise(replies)
            steps = steps + steps_ac
        names = [a.get("callsign") or a["icao"] for a in aircraft]
        uav = [a.get("callsign") or a["icao"] for a in aircraft if a.get("uav")]
        res = result("success", "Mode S/ADS-B receiver (this package) + pyModeS", steps, f,
                     label="Mode S / ADS-B 1090 MHz", confidence=0.97,
                     modulation="PPM 1 Mb/s (pulse position)",
                     evidence=f"{len(msgs)} Mode S frame(s) with valid parity, aircraft {names}"
                              + (f"; UAV (emitter category B6): {uav}" if uav else "")
                              + (f"; {len(replies)} Mode A/C repl(ies)" if replies else ""),
                     span=(F_ADSB - 1e6, F_ADSB + 1e6))
        return [receiver_signal([(m["t_s"], m["t_s"] + (8 + m["bits"]) * 1e-6) for m in msgs]
                                + [(r["t_s"], r["t_s"] + modeac.SPI_US * 1e-6) for r in replies],
                                F_ADSB - 1e6, F_ADSB + 1e6, res, None, self.name)]

    def _mode_ac_signal(self, replies, steps):
        """Mode A/C replies without any Mode S frame: their own signal (research)."""
        f = modeac.summarise(replies)
        rep = [g for g in f["replies"] if g["repeated"]]
        special = sorted({f"{g['code']} {g['special_code']}" for g in f["replies"]
                          if g.get("special_code") and g["repeated"]})
        codes = [g["code"] for g in (rep or f["replies"])[:6]]
        label = ("SSR Mode A/C replies" if rep else "SSR Mode A/C replies-like") + \
            (f" (squawk {', '.join(special)})" if special else "")
        res = result("partial", "Mode A/C reply decoder (this package; research)", steps, {"mode_ac": f},
                     "Mode A/C replies carry no parity: codes are verified only by repetition",
                     label=label, confidence=0.8 if rep else 0.6, modulation="pulse code (SSR Mode A/C)",
                     evidence=f"{f['n_replies']} reply(ies) on the 1.45 us grid between framing pulses 20.3 us "
                              f"apart, {f['n_distinct']} distinct, {f['n_repeated']} repeated; codes {codes}",
                     span=(F_ADSB - 1e6, F_ADSB + 1e6))
        return receiver_signal([(r["t_s"], r["t_s"] + modeac.SPI_US * 1e-6) for r in replies],
                               F_ADSB - 1e6, F_ADSB + 1e6, res, None, self.name)



class UatDecoder(Decoder):
    """UAT 978 MHz ADS-B downlink (research, PLAN G14): capture-level, like Mode S."""
    name = "uat"
    protocol = "UAT 978 MHz ADS-B"
    capture_level = True

    def can_decode(self, sig, ctx):
        from .uat import F_UAT
        return sig.f0 - 1e6 <= F_UAT <= sig.f1 + 1e6

    def capture_signals(self, ctx):
        from . import uat as U
        cap = ctx.cap
        if not (cap.fc - 0.45 * cap.fs < U.F_UAT < cap.fc + 0.45 * cap.fs):
            return []
        y, _ = ctx.ch.slice(U.F_UAT, 0.8e6, U.FS)
        msgs = U.demod(y)
        if not msgs:
            ctx.cache["uat_note"] = "no UAT downlink message passed Reed-Solomon in this capture"
            return []
        steps = ["channelise 978 MHz +-0.8 MHz at 4 samples / bit (4.1667 MSps)",
                 "phase change per bit, 36-bit sync 0xEACDDA4E2 (<= 4 bit errors) at the best timing phase, slicer "
                 "centred on the sync's levels (CFO)",
                 "Reed-Solomon GF(256) 0x187, fcr 120: Long RS(48,34) when the payload type is not 0, else Basic "
                 "RS(30,18) (as FlightAware dump978)",
                 "DO-282B fields: header, state vector, Mode Status (emitter category, callsign)"]
        ac = {}
        for m in msgs:
            a = ac.setdefault(m["address"], {"address": m["address"], "qualifier": m["address_qualifier"], "n": 0})
            a["n"] += 1
            for k in ("lat", "lon", "alt_ft", "callsign", "emitter_category", "uav"):
                if m.get(k) is not None:
                    a[k] = m[k]
        aircraft = list(ac.values())
        uav = [a.get("callsign") or a["address"] for a in aircraft if a.get("uav")]
        f = {"n_messages": len(msgs), "aircraft": aircraft, "messages": msgs}
        res = result("success", "UAT 978 receiver (this package; research)", steps, f,
                     label="UAT 978 MHz ADS-B" + (", UAV" if uav else ""), confidence=0.95,
                     modulation="CPFSK 1.041667 Mb/s, h 0.6",
                     evidence=f"{len(msgs)} RS-corrected UAT message(s), aircraft "
                              f"{[a.get('callsign') or a['address'] for a in aircraft]}"
                              + (f"; UAV (emitter category 14): {uav}" if uav else ""),
                     span=(U.F_UAT - 0.8e6, U.F_UAT + 0.8e6))
        return [receiver_signal([(m["t_s"], m["t_s"] + (36 + (8 * 48 if m["frame"] == "long" else 8 * 30)) / U.RATE)
                                 for m in msgs], U.F_UAT - 0.8e6, U.F_UAT + 0.8e6, res, None, self.name)]


PLUGINS = [AdsbDecoder(), UatDecoder()]
