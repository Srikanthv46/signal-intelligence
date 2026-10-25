"""Pipeline plugin: Herelink and unknown-OFDM digital video links (wraps .detector)."""
from __future__ import annotations

import re
from collections import Counter

from spectrum_intel.core.bands_india import DIGITAL_LINK_RANGES, band_keys
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean, _in_ranges
from spectrum_intel.decode.base import result

CELLULAR = re.compile(r"^[bn]\d")    # LTE / NR band keys in core.bands_india
HERELINK_MIN_HZ = 2.400e9          # below: LTE band 40 (same PHY), left to the LTE decoder
NUMEROLOGY_NAMES = {"lte_like": "LTE-like: DJI OcuSync / O3, CubePilot, Autel", "dji_like": "DJI OcuSync-like",
                    "autel_skylink_like": "Autel SkyLink-like", "herelink": "Herelink-like",
                    "dji_o4_like": "DJI O4-like", "lightbridge_like": "DJI Lightbridge / WiMAX-like",
                    "walksnail_like": "Walksnail / WiMAX-like", "pmp450_like": "Cambium PMP 450-like",
                    "hdzero_like": "HDZero-like"}


class DigitalVtxDecoder(EmissionReceiver):
    name = "digital_vtx"
    protocol = "Herelink / HDZero / unknown OFDM video link"
    key = "digital_vtx"
    method = "receivers.digital_vtx (research receiver)"
    steps = ["Herelink: LTE PSS (ZC 25/29/34) fold on the 5 ms half-frame over the band, carrier search, "
             "SSS decode -> PCI, CP-phase CFO, subframe activity -> control / video / uplink",
             "HDZero (27 MHz mode): on each HDZero channel in the capture, 1024-FFT at 27 MSps slice, CP-fold symbol "
             "timing and fractional CFO, continual-pilot coherence vs the data carriers beside them over the integer "
             "carrier offsets -> channel, ppm, per-frame preamble -> video frame rate",
             "any other CP-OFDM link (Tu >= 10 us): coherent lag-Tu autocorrelation scan (time-AGC'd), symbol-rate "
             "line / per-window fold confirmation, CP cross-spectrum band, bursts on local symbol grids, TDD frame "
             "completion, band and burst extent from in-band power seeded by the CP evidence; numerology table -> "
             "drone family / infrastructure / unknown",
             "hopping narrowband CP-OFDM (RC uplinks): windowed CP-fold scan at the capture's and the drone families' "
             "numerologies, envelope edges, on-vs-off band fit, CP-fold verification, frame-cadence tracking"]

    def run_detector(self, ctx):
        from .detector import detect
        return detect(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc, {"report_wifi": False})

    def in_range(self, lo, hi):
        return _in_ranges(lo, hi, DIGITAL_LINK_RANGES)

    def keep(self, e):
        if e.family == "herelink":
            return e.center_mhz * 1e6 >= HERELINK_MIN_HZ and e.confidence >= 0.5
        if e.family == "hdzero":
            return True
        if e.family != "digital_vtx":
            return False
        # inside a licensed cellular band (e.g. LTE B40 TDD at 2300-2400 MHz) 15/30 kHz OFDM is
        # the operator's carrier: leave it to the lte/nr receivers and the classifier
        if any(CELLULAR.match(k) for k in band_keys(e.center_mhz * 1e6, e.center_mhz * 1e6)):
            return False
        # 802.11-family numerology at a scaled clock (12.8 us = 802.11a at 5 MHz) is the Wi-Fi receivers'
        return e.features.get("role") != "wifi"

    def group_key(self, e):
        if e.family in ("herelink", "hdzero"):
            return (e.family, "", round(e.center_mhz, 0))
        if e.label.startswith("ofdm_hop"):
            return (e.family, e.label, e.features.get("channel_mhz"))      # one signal per hop channel
        return (e.family, e.label, e.features.get("carrier_mhz"))

    def arbitrate(self, signals, ctx):
        """A blind OFDM link that a protocol receiver identified -- a DJI video link or Herelink of
        the same subcarrier spacing over more than half of its band -- is that receiver's signal:
        it is dropped and recorded on the host (explained_receiver_detections)."""
        gen = [s for s in signals if s.features.get("decoder") == self.name
               and not s.label.startswith(("Herelink", "HDZero"))]
        if not gen:
            return signals
        hosts = []
        for s in signals:
            if s.features.get("decoder") == "dji":
                scs = {e.features.get("scs_khz") for e in s.features.get("_emissions", [])
                       if e.label != "uplink"} - {None}
            elif s.features.get("decoder") == self.name and s.label.startswith("Herelink"):
                scs = {15.0}
            elif s.features.get("decoder") == self.name and s.label.startswith("HDZero"):
                scs = {26.367}
            else:
                continue
            hosts.append((s, scs))
        # a blind hop that a DJI RC uplink hop covers (same time, same band) is the DJI receiver's; hops
        # that mostly sit inside a DJI video link's bursts and band are slices of that video
        ups = [s for s in signals if s.features.get("decoder") == "dji" and "uplink" in s.label]
        vids = [s for s in signals if s.features.get("decoder") == "dji" and "video" in s.label]
        keep = []
        for s in signals:
            if s in gen and s.decoded.get("fields", {}).get("hop_us") is not None:
                host = next((u for u in ups if min(u.f1, s.f1) - max(u.f0, s.f0) > 0.5 * min(u.bw, s.bw)
                             and any(min(a1, b1) - max(a0, b0) > 0 for a0, a1 in _times(u) for b0, b1 in _times(s))),
                            None)
                note = "the same hops found blind (numerology, band, timing): corroborates it"
                if host is None:
                    host = next((v for v in vids if min(v.f1, s.f1) - max(v.f0, s.f0) > 0.5 * s.bw and 2 * sum(
                        any(min(b1, a1) - max(b0, a0) > 0.5 * (b1 - b0) for a0, a1 in _times(v)) for b0, b1 in _times(s))
                        > len(_times(s))), None)
                    note = "narrowband slices of this link's bursts, found blind: part of it"
                if host is not None:
                    host.detection.setdefault("explained_receiver_detections", []).append(
                        {"decoder": self.name, "label": s.label, "f_low_hz": round(s.f0), "f_high_hz": round(s.f1),
                         "note": note})
                    continue
                keep.append(s)
                continue
            if s in gen:
                F = s.decoded.get("fields", {})
                host = next((h for h, scs in hosts if min(h.f1, s.f1) - max(h.f0, s.f0) > 0.5 * min(h.bw, s.bw)
                             and any(abs(F.get("scs_khz", 0) / v - 1) < 0.03 for v in scs)), None)
                if host is not None:
                    host.detection.setdefault("explained_receiver_detections", []).append(
                        {"decoder": self.name, "label": s.label, "f_low_hz": round(s.f0), "f_high_hz": round(s.f1),
                         "note": "the same OFDM link found blind (numerology, band, bursts): corroborates it"})
                    continue
            keep.append(s)
        return keep

    def summarise(self, em, sig):
        e = max(em, key=lambda x: x.confidence)
        F = e.features
        labels = Counter(x.label for x in em)
        keys = ("pci", "nid1", "carrier_mhz", "cfo_hz", "cfo_ppm", "bw_mhz", "scs_khz", "tu_us", "cp_us",
                "symbol_us", "numerology", "frame_period_ms", "burst_us", "clock_ppm", "pss_instances", "role",
                "matches", "cp_fraction", "duty_cycle", "symbol_grids", "line_z", "strong_windows", "link_snr_db",
                "channel_mhz", "hop_bw_mhz", "hop_us", "n_hops", "n_channels", "hop_period_ms", "paired_link_mhz",
                "tx_antennas", "system_candidates", "channel", "bandwidth_mode", "fft_size", "sample_rate_msps",
                "used_subcarriers", "continual_pilots", "pilot_coherence", "data_coherence", "cp_fold_contrast",
                "n_symbols", "preamble_ms", "frame_symbols", "frame_rate_fps", "band_context")
        fields = _clean({k: F.get(k) for k in keys if F.get(k) is not None})
        fields.update(n_bursts=len(em), labels=dict(labels))
        if e.label.startswith("ofdm_hop"):            # model level from the hop timing (research, PLAN M6)
            from .catalog import hop_model
            mb = hop_model(F.get("hop_bw_mhz"), F.get("hop_us"), F.get("hop_period_ms"))
            if mb:
                fields["model"] = mb
        if F.get("system_candidates"):                # model level (research, PLAN M3 / M4)
            from .catalog import model_candidates
            mb = model_candidates(F["system_candidates"], F.get("carrier_mhz"), F.get("bw_mhz"),
                                  (F.get("tx_antennas") or {}).get("tx_antennas"))
            if mb:
                fields["model"] = mb
        if e.family == "hdzero":
            from .catalog import hdzero_model
            fields["model"] = hdzero_model(F.get("bandwidth_mode"))
            status = "partial"
            reason = ("HDZero PHY synchronised (symbol timing, carrier, pilots, video frame); the video samples are "
                      "mapped straight onto the subcarriers (JSCC) with no published format: not decoded")
            label = "HDZero digital FPV video link"
            mod = (f"OFDM {F.get('scs_khz')} kHz SCS (1024 FFT at 27 MSps), CP {F.get('cp_us')} us, "
                   f"{F.get('used_subcarriers')} subcarriers")
            fps = F.get("frame_rate_fps")
            ev = (f"HDZero continual pilots coherent on channel {F.get('channel')} ({F.get('carrier_mhz')} MHz, "
                  f"{F.get('cfo_ppm')} ppm): pilot coherence {F.get('pilot_coherence')} vs {F.get('data_coherence')} "
                  f"on the data carriers beside them over {F.get('n_symbols')} symbols, CP fold contrast "
                  f"{F.get('cp_fold_contrast')}, {F.get('bw_mhz')} MHz, SNR {F.get('link_snr_db')} dB"
                  + (f"; per-frame preamble every {F.get('frame_period_ms')} ms ({fps} fps)" if fps else
                     f"; preamble at {F.get('preamble_ms')} ms" if F.get("preamble_ms") else ""))
        elif e.family == "herelink":
            pci = F.get("pci")
            status = "partial"
            reason = ("LTE-like PHY synchronised (PSS" + ("+SSS -> PCI" if pci is not None else "") +
                      "); video payload is proprietary/encrypted")
            label = "Herelink (CubePilot) drone video/control link"
            mod = "OFDM 15 kHz SCS (LTE Rel-8 DL PHY)"
            cm = F.get("carrier_mhz") or e.center_mhz
            ism = "2.4 GHz" if cm < 3000 else "5.8 GHz" if cm >= 5650 else "5 GHz"
            ev = (f"LTE PSS at {F.get('carrier_mhz')} MHz in the {ism} band with Herelink TDD subframe "
                  f"rhythm ({dict(labels)})" + (f", PCI {pci}" if pci is not None else ""))
        elif e.label.startswith("ofdm_hop"):
            status = "partial"
            scs = F.get("scs_khz", 0)
            sc = f"{scs:.2f}".rstrip("0").rstrip(".")
            num, role = F.get("numerology"), F.get("role")
            hopping = (F.get("n_channels") or 0) >= 2
            kind = "hopping" if hopping else "narrowband burst"
            if role == "drone":
                label = f"OFDM {'hopping ' if hopping else ''}RC uplink-like ({NUMEROLOGY_NAMES.get(num, num)} numerology, {sc} kHz SCS)"
            elif role == "infrastructure":
                label = f"OFDM fixed-wireless uplink-like ({NUMEROLOGY_NAMES.get(num, num)} numerology, {sc} kHz SCS)"
            else:
                label = f"unknown OFDM {kind} link ({sc} kHz SCS)"
            reason = (f"narrowband CP-OFDM bursts ({'hopping between channels' if hopping else 'on one channel'}), "
                      "found and measured blindly; no protocol receiver identified them")
            mod = f"OFDM {sc} kHz SCS, CP {F.get('cp_us')} us, ~{F.get('hop_bw_mhz')} MHz hops"
            pl = F.get("paired_link_mhz") or []
            ev = (f"blind CP-OFDM hops: Tu {F.get('tu_us')} us, CP {F.get('cp_us')} us, {len(em)} hop(s) on this "
                  f"channel ({F.get('channel_mhz')} MHz), {F.get('n_hops')} hops on {F.get('n_channels')} channel(s) in "
                  f"the capture, {F.get('hop_us')} us x {F.get('hop_bw_mhz')} MHz each, hop period "
                  f"{F.get('hop_period_ms')} ms" + (f"; same numerology as the OFDM link at {pl} MHz" if pl else ""))
        else:
            status = "partial"
            scs = F.get("scs_khz", 0)
            sc = f"{scs:.2f}".rstrip("0").rstrip(".")
            role, num = F.get("role"), F.get("numerology")
            if role == "drone":
                label = f"OFDM drone video link-like ({NUMEROLOGY_NAMES.get(num, num)} numerology, {sc} kHz SCS)"
            elif role == "infrastructure":
                label = f"OFDM fixed-wireless link-like ({NUMEROLOGY_NAMES.get(num, num)} numerology, {sc} kHz SCS)"
            else:
                label = f"unknown OFDM {'TDD ' if F.get('frame_period_ms') else ''}link ({sc} kHz SCS)"
            if F.get("band_context"):                       # outside the ISM / FPV bands (PLAN G8)
                label += " in the " + F["band_context"].split(" (")[0]
            reason = ("CP-OFDM link found and measured blindly; the numerology resembles a known system but no "
                      "protocol receiver identified it" if role in ("drone", "infrastructure") else
                      "CP-OFDM link found and measured blindly; numerology unknown, no protocol to decode")
            mod = f"OFDM {sc} kHz SCS, CP {F.get('cp_us')} us"
            ev = (f"blind CP-OFDM: Tu {F.get('tu_us')} us, CP {F.get('cp_us')} us ({F.get('cp_fraction')} Tu), "
                  f"{F.get('bw_mhz')} MHz at {F.get('carrier_mhz')} MHz, {len(em)} burst(s), duty {F.get('duty_cycle')}, "
                  f"frame period {F.get('frame_period_ms')} ms, symbol-rate line z {F.get('line_z')}, "
                  f"{F.get('strong_windows')} strong 1 ms windows, SNR {F.get('link_snr_db')} dB")
            if F.get("band_context"):
                ev += f"; band: {F['band_context']}"
            tx = F.get("tx_antennas") or {}
            if tx.get("tx_antennas") == 2:                  # inferred from the data statistics: "-like"
                scheme = "STBC" if tx["verdict"] == "2tx_stbc_like" else "SFBC"
                label = label[:-1] + f", 2-TX {scheme}-like)" if label.endswith(")") else label + f" (2-TX {scheme}-like)"
                mod += f", 2 transmit antennas ({scheme}-like)"
            if tx.get("verdict"):
                ev += f"; transmit streams: {tx['verdict']}" + (f" ({tx['note']})" if tx.get("note") else "")
            sc = F.get("system_candidates") or []
            if sc:
                ev += (f"; most like {sc[0]['system']} (research catalog, matched {', '.join(sc[0]['matched'])}"
                       + (f"; conflicts {', '.join(sc[0]['conflicts'])}" if sc[0]['conflicts'] else "") + ")")
        return result(status, self.method, self.steps, fields, reason, label=label,
                      confidence=e.confidence, modulation=mod, evidence=ev,
                      span=(min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6))


def _times(s):
    return [(e.start_us, e.end_us) for e in s.features.get("_emissions", [])]


PLUGINS = [DigitalVtxDecoder()]
