"""Pipeline plugin: 802.11b DSSS/CCK PPDUs (Wi-Fi OFDM is in plugin_ofdm.py).

One signal per 2.4 GHz channel with decoded PPDUs. Its fields hold the PLCP statistics, a table of the
BSSs (from beacons / probe responses and the frames that name them), the stations heard, and one line per
frame. Barker-coherent energy without a decodable PLCP header is reported as such (partial), not as unknown.
"""
from __future__ import annotations

from collections import Counter

import numpy as np

from spectrum_intel.core.bands_india import WIFI_24
from spectrum_intel.decode.base import Decoder, overlaps, receiver_signal, result

from .dsss import scan
from .plugin_ofdm import drone_label

STEPS = ["channelise each 2.4 GHz channel: 22 MSps (detection) and 44 MSps (4 samples/chip, decoding)",
         "detect: |Barker-11 matched filter|^2 summed over 16 symbols at each of the 22 sample phases; DSSS where "
         "the best phase stands out (peak/mean > 2.8, noise < 2.6); runs split where the best phase moves",
         "acquire SFD-first: Barker symbols on a 1 us grid at every sample phase, DBPSK, descramble x^7+x^4+1, "
         "SFD (long 0xF3A0 / short 0x05CF, <= 1 bit error) after a SYNC run; sidelobe hits merged",
         "track each PPDU from its own SYNC (early/late loop gated on its amplitude); rebuild the scrambled "
         "SYNC+SFD exactly (best of 128 scrambler states) -> data-aided CFO, phase, SNR",
         "PLCP header coherent (phase-tracked), differential fallback; CRC-16",
         "PSDU: 1/2 Mb/s coherent Barker symbols; 5.5/11 Mb/s CCK on a CFO-locked chip grid through a 9-tap "
         "LS equaliser trained on preamble + header chips, ML codeword + phi1 decisions",
         "descramble, CRC-32 FCS; on failure FCS-guided correction of <= 4 of the 20 least reliable symbol "
         "decisions, then STFT interference excision and a re-decode",
         "MAC dissection (mac80211.py): addresses and roles, IEs, security, capabilities, WPS / P2P "
         "identity, EAPOL / ARP / IP headers when unencrypted; vendor from the IEEE registry",
         "dedup: a PPDU seen on neighbouring channels is kept where its Barker amplitude is highest; "
         "a beacon's DS parameter names the channel"]


class WifiDsssDecoder(Decoder):
    name = "wifi_dsss"
    protocol = "IEEE 802.11b DSSS/CCK"
    capture_level = True

    def can_decode(self, sig, ctx):
        return 2401e6 <= sig.fc <= 2495e6 and sig.bw >= 3e6

    def capture_signals(self, ctx):
        if not (ctx.cap.fc - ctx.cap.fs / 2 < 2495e6 and ctx.cap.fc + ctx.cap.fs / 2 > 2401e6):
            return []
        pk = [p for p in scan(ctx) if p.get("plcp_crc_ok")]
        groups = {}
        for p in pk:
            groups.setdefault(p["channel"], []).append(p)
        out = []
        for ch, g in groups.items():
            res = self.summarise(g, len(pk), ctx)
            f0, f1 = res["span"]
            # signal SNR comes from the energy detections this signal absorbs (same scale as every other
            # signal); the receiver's in-band SNR per PPDU is in fields (snr_db_inband_median, frames)
            out.append(receiver_signal([(p["start_s"], p["end_s"]) for p in g], f0, f1, res, None, self.name))
        return out

    def decode(self, sig, ctx):
        pk = scan(ctx)
        mine = [p for p in pk if p.get("plcp_crc_ok") and
                overlaps(sig, p["start_s"], p["end_s"], p["ch_freq_mhz"] * 1e6 - 11e6, p["ch_freq_mhz"] * 1e6 + 11e6)]
        if mine:
            return self.summarise(mine, len(pk), ctx)
        regs = [r for r in ctx.cache.get("wifi_dsss_unexplained", [])
                if overlaps(sig, r["t0"], r["t1"], r["ch_freq_mhz"] * 1e6 - 11e6, r["ch_freq_mhz"] * 1e6 + 11e6)]
        if regs:
            return self.barker_only(regs, ctx)
        return result("not_possible", "802.11b DSSS receiver (this package)", STEPS,
                      reason="no Barker-11 spreading (coherence < 0.16) and no PLCP header in this signal's span")

    def barker_only(self, regs, ctx):
        r = max(regs, key=lambda r: r["coherence"])
        early = min(x["t0"] for x in regs) < 60e-6
        dur = sum(x["t1"] - x["t0"] for x in regs)
        fields = {"channel": r["channel"], "center_mhz": r["ch_freq_mhz"], "n_regions": len(regs),
                  "barker_coherence": r["coherence"], "duration_us": round(dur * 1e6),
                  "regions": [{"t_start_s": round(x["t0"], 6), "t_end_s": round(x["t1"], 6), "channel": x["channel"],
                               "coherence": x["coherence"]} for x in regs[:50]]}
        if early:
            lab, conf = "Wi-Fi 802.11b DSSS-like (preamble before capture start)", 0.7
            why = "Barker-11 spreading at 1 Msym/s from the start of the capture: the PPDU began before it"
        else:
            lab, conf = "Barker-11 DSSS (no 802.11b PLCP header)-like", 0.55
            why = ("Barker-11 spreading at 1 Msym/s but no SFD / valid PLCP header: the preamble was lost "
                   "(collision, weak) or the framing is not 802.11b")
        f = r["ch_freq_mhz"] * 1e6
        return result("partial", "802.11b DSSS receiver (this package)", STEPS, fields, reason=why, label=lab,
                      confidence=conf, modulation="DSSS Barker-11 1 Msym/s", span=(f - 11e6, f + 11e6),
                      evidence=f"Barker-11 matched-filter coherence {r['coherence']:.2f} (noise ~0.05) over "
                               f"{dur * 1e6:.0f} us on ch{r['channel']}; {why}")

    def summarise(self, mine, n_all, ctx=None):
        steps = STEPS + [f"{n_all} PLCP headers with valid CRC-16 in capture; {len(mine)} attributed to this signal"]
        if not mine:
            return result("not_possible", "802.11b DSSS receiver (this package)", steps,
                          reason="no Barker-spread PLCP header with valid CRC-16 in this signal's span")
        fcs = [p for p in mine if p.get("fcs_ok")]
        macs = [p["mac"] for p in fcs]
        ch = Counter(p["channel"] for p in mine).most_common(1)[0][0]
        bss, stations = _bss_table(fcs), _station_table(fcs)
        snr = [p["snr_db"] for p in mine if p.get("snr_db") is not None]
        fields = {"channel": ch, "center_mhz": WIFI_24[ch] / 1e6,
                  "n_ppdus": len(mine), "n_fcs_ok": len(fcs),
                  "n_fcs_corrected": sum(1 for p in fcs if p.get("symbols_corrected")),
                  "rates_mbps": dict(Counter(p["rate_mbps"] for p in mine)),
                  "preambles": dict(Counter(p["preamble"] for p in mine)),
                  "channels": dict(Counter(p["channel"] for p in mine)),
                  "snr_db_inband_median": round(float(np.median(snr)), 1) if snr else None,
                  "cfo_hz_median": int(np.median([p["cfo_hz"] for p in mine])),
                  "bssids": sorted({m.get("bssid") for m in macs if m.get("bssid")}),
                  "ssids": sorted({m["ssid"] for m in macs if m.get("ssid")}),
                  "frame_types": dict(Counter(m.get("subtype") for m in macs)),
                  "n_frames": len(fcs),
                  "bss": bss, "stations": stations,
                  "frames": [_frame_line(p) for p in mine[:300]],
                  "ppdus": mine[:200]}
        rids = [m["remote_id"] for m in macs if m.get("remote_id")]
        if rids:                                     # research, PLAN M2: manufacturer from the Remote ID serial
            from spectrum_intel.core.drone_models import remote_id_model
            fields["model"] = remote_id_model(rids[-1])
        if ctx is not None:                          # research, PLAN W1
            from .plugin_ofdm import CTRL_MAX_BYTES, field_summary
            frames = [{"t0": p["start_s"], "t1": p["end_s"], "dur_us": (p["end_s"] - p["start_s"]) * 1e6,
                       "bytes": p.get("psdu_bytes"), "ta": (p.get("mac") or {}).get("ta") or (p.get("mac") or {}).get("sa"),
                       "macs": [p["mac"]] if p.get("fcs_ok") and p.get("mac") else [],
                       "ctrl": (p.get("mac") or {}).get("type") == "ctrl" or (not p.get("fcs_ok") and
                                                                              (p.get("psdu_bytes") or 99) <= CTRL_MAX_BYTES),
                       "video": bool(((p.get("mac") or {}).get("ipv4") or {}).get("video"))} for p in mine]
            fields.update(field_summary(frames, len(ctx.cap.iq) / ctx.cap.fs))
        names = ", ".join(f"'{b['ssid']}'" for b in bss if b.get("ssid"))[:200]
        return result("success" if fcs else "partial", "802.11b DSSS/CCK receiver (this package)", steps, fields,
                      "" if fcs else "PLCP headers decoded (CRC-16 ok) but no PSDU passed FCS",
                      label=f"Wi-Fi 802.11b DSSS ch{ch}" + drone_label(fields)[0], confidence=0.97 if fcs else 0.9,
                      modulation="DSSS Barker-11 DBPSK/DQPSK, CCK", span=(WIFI_24[ch] - 11e6, WIFI_24[ch] + 11e6),
                      evidence=f"{len(mine)} 802.11b PLCP headers with valid CRC-16"
                               + (f"; {len(fcs)} MPDUs pass FCS" if fcs else "")
                               + (f"; {len(bss)} BSS: {names}" if bss else "") + drone_label(fields)[1])


def _frame_line(p):
    m = p.get("mac") or {}
    return {"t_ms": round(p["start_s"] * 1e3, 3), "rate_mbps": p.get("rate_mbps"), "preamble": p.get("preamble"),
            "length_us": p.get("length_us"), "snr_db": p.get("snr_db"), "cfo_khz": round(p["cfo_hz"] / 1e3, 1),
            "fcs_ok": bool(p.get("fcs_ok")), "symbols_corrected": p.get("symbols_corrected", 0),
            "summary": m.get("summary") or ("FCS failed" if p.get("psdu_bytes") else "header only")}


_BSS_KEYS = ("ssid", "ssid_hidden", "security", "phy", "ap_uptime", "beacon_interval_ms", "capabilities",
             "rates_mbps", "basic_rates_mbps", "country", "bss_load", "ht", "ht_operation", "vht_operation",
             "wps", "p2p", "mobility_domain", "radio_measurement_802_11k", "interworking_802_11u",
             "vendor_elements", "erp", "tim")


def _bss_table(frames):
    """One row per BSSID: what its beacons / probe responses say, and the traffic that names it."""
    rows: dict = {}
    for p in frames:
        m = p["mac"]
        b = m.get("bssid")
        if not b or b == "ff:ff:ff:ff:ff:ff":
            continue
        r = rows.setdefault(b, {"bssid": b, "vendor": m.get("bssid_vendor"), "n_frames": 0, "frame_types": Counter(),
                                "stations": set(), "snr_db": []})
        r["n_frames"] += 1
        r["frame_types"][m.get("subtype")] += 1
        if p.get("snr_db") is not None and m.get("ta") == b:
            r["snr_db"].append(p["snr_db"])
        if m.get("subtype") in ("beacon", "probe-resp"):
            r["channel"] = m.get("ds_channel") or p["channel"]
            r["cfo_hz"] = p["cfo_hz"]
            for k in _BSS_KEYS:
                if k in m:
                    r[k] = m[k]
        for k in ("sa", "da", "ta", "ra"):
            v = m.get(k)
            if v and v != b and not int(v[:2], 16) & 1:
                r["stations"].add(v)
    out = []
    for r in rows.values():
        r["frame_types"] = dict(r["frame_types"])
        r["stations"] = sorted(r["stations"])
        r["snr_db"] = round(float(np.median(r["snr_db"])), 1) if r["snr_db"] else None
        out.append(r)
    return sorted(out, key=lambda r: (r.get("ssid") is None, -r["n_frames"]))


def _station_table(frames):
    """Non-AP devices heard (transmitters that are not a BSSID), with what they did."""
    aps = {p["mac"].get("bssid") for p in frames if p["mac"].get("subtype") in ("beacon", "probe-resp")}
    rows: dict = {}
    for p in frames:
        m = p["mac"]
        tx = m.get("sa") or m.get("ta")
        if not tx or tx in aps or int(tx[:2], 16) & 1:
            continue
        r = rows.setdefault(tx, {"mac": tx, "vendor": m.get("sa_vendor") or m.get("ta_vendor"),
                                 "randomized": bool(int(tx[:2], 16) & 2), "n_frames": 0, "frame_types": Counter(),
                                 "probed_ssids": set(), "bssid": None})
        r["n_frames"] += 1
        r["frame_types"][m.get("subtype")] += 1
        if m.get("subtype") == "probe-req" and m.get("ssid") not in (None, "<broadcast>"):
            r["probed_ssids"].add(m["ssid"])
        if m.get("bssid") and m.get("bssid") != "ff:ff:ff:ff:ff:ff" and m.get("type") == "data":
            r["bssid"] = m["bssid"]
        name = (m.get("p2p") or {}).get("device_name") or (m.get("wps") or {}).get("device_name")
        if name:
            r["device_name"] = name
        if m.get("phy"):
            r["phy"] = m["phy"]
    out = []
    for r in rows.values():
        r["frame_types"] = dict(r["frame_types"])
        r["probed_ssids"] = sorted(r["probed_ssids"])
        out.append(r)
    return sorted(out, key=lambda r: -r["n_frames"])


PLUGINS = [WifiDsssDecoder()]
