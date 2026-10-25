"""Pipeline plugin: 315 / 433 MHz remotes and TPMS (wraps .detector; research, signal gaps S6)."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean

from .detector import BANDS


def _overlap(lo, hi):
    return any(lo < b and hi > a for a, b in BANDS)


MOD = {"Continental / VDO": "FSK, Manchester, 52 us", "Toyota / Pacific PMV-C210": "FSK, differential Manchester, 52 us"}


class IsmDevicesDecoder(EmissionReceiver):
    name = "ism_devices"
    protocol = "315 / 433 MHz remotes and TPMS"
    key = "ism_devices"
    method = "receivers.ism_devices: OOK pulse timing / FSK symbols -> EV1527, PT2262, KeeLoq, TPMS (research)"
    steps = ["candidates: 314-316 and 433.05-434.79 MHz at >= 2 MSps, 256-point frames; per frame the bin furthest "
             "over its median across frames (>= 12 dB); frequencies grouped at 100 kHz gaps -> up to 4 carriers",
             "per carrier +-100 kHz at 1 MSps: envelope (8-sample average) on above the geometric mean of its 10th / "
             "99.5th percentiles (>= 12 dB apart), runs < 30 us debounced -> pulses (mark, space)",
             "OOK: EV1527 / PT2262 (sync 1 + 31 Te, 24 PWM bits of 4 Te, repeats compared; PT2262 tri-state when "
             "every pair is 00 / 11 / 01); KeeLoq HCS (12 warm-up pulses, 10 Te header, 66 PWM bits of 3 Te: clear "
             "serial, buttons, battery, repeat; the encrypted hopping code is not reported)",
             "FSK: a mark >= 2 ms with a two-level discriminator (deviation >= 10 kHz) -> 52 us symbols -> TPMS "
             "Continental / VDO (Manchester, CRC-8 0x07 / 0xAA) and Toyota PMV-C210 (differential Manchester, "
             "CRC-8 0x07 / 0x80 + inverted pressure), both polarities",
             "nothing decoded: >= 16 OOK pulses of 50-2000 us or an FSK burst -> an unrecognised device ('-like')"]

    def run_detector(self, ctx):
        from .detector import detect
        return detect(ctx.ch)

    def in_range(self, lo, hi):
        return _overlap(lo, hi)

    def can_decode(self, sig, ctx):
        return _overlap(sig.f0, sig.f1)

    def arbitrate(self, signals, ctx):
        """A device this receiver decoded (TPMS CRC, repeated fixed code, KeeLoq) explains an unverified receiver
        guess (no success decode, e.g. an FSK rate 'SiK-like') whose centre lies in its span and whose bursts
        overlap its frames in time: dropped and recorded on the host (explained_receiver_detections)."""
        hosts = [s for s in signals if s.features.get("decoder") == self.name and s.confidence >= 0.75]
        if not hosts:
            return signals
        keep = []
        for s in signals:
            host = None
            if s.features.get("decoder") not in (None, self.name) and s.decoded.get("status") != "success":
                host = next((h for h in hosts if h.f0 <= s.fc <= h.f1 and any(
                    min(a.t1, b.t1) > max(a.t0, b.t0) for a in s.bursts for b in h.bursts)), None)
            if host is None:
                keep.append(s)
                continue
            host.detection.setdefault("explained_receiver_detections", []).append(
                {"decoder": s.features.get("decoder"), "label": s.label, "f_low_hz": round(s.f0),
                 "f_high_hz": round(s.f1), "note": "the same burst as this decoded device: attributed to it"})
        return keep

    def group_key(self, e):
        d = e.decoded or {}
        return (e.label, d.get("id") or d.get("code") or d.get("serial"), round(e.center_mhz, 2))

    def summarise(self, em, sig):
        e = max(em, key=lambda e: e.confidence)
        d = dict(e.decoded or {})
        F = e.features
        span = (min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6)
        base = {"carrier_mhz": F["carrier_mhz"], "envelope_db": F["envelope_db"]}
        if e.label == "tpms":
            fields = _clean({**d, "frames": len(em), **base})
            return result("success", self.method, self.steps, fields, "",
                          label=f"TPMS ({d['protocol']}, ID {d['id']})", confidence=e.confidence,
                          modulation=MOD[d["protocol"]],
                          evidence=f"{len(em)} CRC-8 verified frame(s); {d['pressure_kpa']} kPa, "
                                   f"{d['temperature_c']} C", span=span)
        if e.label == "ook_fixed_code":
            kind = "PT2262 tri-state " + d["tristate"] if "tristate" in d else "EV1527"
            return result("partial", self.method, self.steps, _clean({**d, **base}),
                          "fixed-code remotes carry no checksum: the repeats are the evidence",
                          label=f"OOK fixed-code remote (EV1527 / PT2262-like, ID {d['id']})",
                          confidence=e.confidence, modulation="OOK PWM",
                          evidence=f"{d['repeats']} identical 24-bit word(s) of {d['words']}, Te {d['te_us']} us; "
                                   f"fits {kind}", span=span)
        if e.label == "keeloq":
            return result("partial", self.method, self.steps, _clean({**d, **base}),
                          "code-hopping remote: the encrypted 32-bit part is neither decoded nor reported",
                          label=f"KeeLoq code-hopping remote-like (serial {d['serial']})", confidence=e.confidence,
                          modulation="OOK PWM",
                          evidence=f"{d['packets']} packet(s) of 12 warm-up pulses + 66 PWM bits, Te {d['te_us']} us; "
                                   f"{d['matching_packets']} with the same serial and button", span=span)
        kind = "OOK" if e.label == "ook_device" else "FSK"
        return result("partial", self.method, self.steps, _clean(dict(F)),
                      "no known protocol matched",
                      label=f"{kind} device-like (unrecognised protocol)", confidence=e.confidence,
                      modulation=kind,
                      evidence=(f"{F['n_pulses']} OOK pulses, marks {F['mark_short_us']}-{F['mark_long_us']} us"
                                if kind == "OOK" else f"{F['fsk_bursts']} FSK burst(s) >= 2 ms"), span=span)


PLUGINS = [IsmDevicesDecoder()]
