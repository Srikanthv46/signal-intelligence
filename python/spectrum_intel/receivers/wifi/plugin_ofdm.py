"""Pipeline plugins: 802.11a/g/n/ac/ax (OFDM) — wrap receivers.wifi.ofdm_rx / ofdm_scaled / payload / links.

wifi_ofdm         Capture-wide scan of every standard 20 MHz channel (L-STF autocorrelation -> L-LTF
                  channel estimate -> L-SIG Viterbi + parity, HT/VHT-SIG CRC); non-HT and single-stream
                  HT (BCC or LDPC, with or without STBC) PPDUs are demodulated to MPDUs and FCS-checked;
                  MAC headers give BSSID/SA/DA and beacons the SSID.
wifi_ofdm_scaled  The same receiver on links that run 802.11 at another clock (5/10 MHz channels, vendor
                  clocks such as 11 MHz), at any frequency, with the spectrum possibly mirrored: the OFDM
                  numerology is measured, then the capture is time-scaled into a standard channel.
Both attribute every PPDU to a transmitter (links.attribute): MAC header, control responses, scrambler
counters, then received power as a weak hint.
"""
from __future__ import annotations

from collections import Counter

import numpy as np

from spectrum_intel.decode.base import Decoder, DecodeContext, overlaps, receiver_signal, result
from spectrum_intel.core.bands_india import band_keys

WIFI_BANDS = {"ism_24", "unii_12", "unii_2e", "ism_58", "its_59", "wifi6e"}


def _records(pkts, cap, clock_hz=20e6):
    """ofdm_rx Packets (real units) -> PPDU records with decoded payloads, seeds and power."""
    from . import payload as WP
    scale = clock_hz / 20e6                                    # real rate = nominal rate * scale
    out = []
    for p in pkts:
        rec = {"t0": p.start_sample / cap.fs, "t1": p.end_sample / cap.fs,
               "f0": p.f_lo_mhz * 1e6, "f1": p.f_hi_mhz * 1e6, "channel": p.channel,
               "center_mhz": p.center_mhz, "fmt": p.fmt, "bw_mhz": p.bandwidth_mhz,
               "rate_mbps": round(p.lsig_rate_mbps * scale, 3), "length": p.lsig_length, "mcs": p.mcs,
               "sig_crc_ok": p.sig_crc_ok, "snr_db": p.snr_db, "cfo_khz": p.cfo_khz,
               "start_us": p.start_us, "dur_us": p.duration_us}
        rx = p.extra.get("_rx")
        if rx and rx[0] is not None:
            st = rx[0]
            seg = st["r"][st["start"]:st["start"] + 320]          # L-STF + L-LTF
            if len(seg):
                rec["power_db"] = round(float(10 * np.log10(np.mean(np.abs(seg) ** 2) + 1e-30)), 1)
        if p.fmt == "NON_HT" and rx and rx[0] is not None:
            try:
                mpdu, ok, seed, info, nsym = WP.decode_nonht_payload(rx[0], p.lsig_rate_mbps, p.lsig_length,
                                                                     rx[1], rx[2])
                rec["fcs_ok"] = bool(ok)
                rec["scrambler_seed"] = int(sum(int(b) << i for i, b in enumerate(seed)))
                rec["seed_ok"] = bool(ok or info.get("_service_ok"))
                if ok:
                    rec["mac"] = {k: v for k, v in info.items() if not k.startswith("_")}
                    rec["macs"] = [rec["mac"]]
            except Exception as e:                      # noqa: BLE001 — record, don't fail the run
                rec["payload_error"] = str(e)[:120]
        elif p.fmt == "HT" and rx and rx[0] is not None:
            e = p.extra
            rec["psdu_bytes"] = e.get("ht_length")
            rec["ht_rate_mbps"] = (round(WP.ht_rate_mbps(p.mcs, bool(e.get("short_gi"))) * scale, 3)
                                   if p.mcs in WP.HT_MCS else None)
            rec["short_gi"] = bool(e.get("short_gi"))
            rec["aggregated"] = bool(e.get("aggregation"))
            rec["coding"] = "LDPC" if e.get("ldpc") else "BCC"
            rec["stbc"] = int(e.get("stbc", 0))
            ok, why = WP.ht_decodable(p.mcs, e.get("ldpc", 0), e.get("stbc", 0), e.get("n_ess", 0),
                                      20 if clock_hz != 20e6 else p.bandwidth_mhz)
            if not ok:
                rec["payload_note"] = why               # say what cannot be decoded, and why
            else:
                try:
                    r = WP.decode_ht_payload(rx[0], p.mcs, e["ht_length"], bool(e.get("short_gi")),
                                             bool(e.get("aggregation")), n_ess=e.get("n_ess", 0),
                                             chan_offset_hz=rx[1], usable_half_hz=rx[2], stbc=e.get("stbc", 0),
                                             ldpc=e.get("ldpc", 0), smoothing=e.get("smoothing", 1))
                    rec["fcs_ok"] = r["n_fcs_ok"] > 0
                    rec["n_subframes"], rec["n_fcs_ok"] = r["n_subframes"], r["n_fcs_ok"]
                    rec["scrambler_seed"] = r["scrambler_seed"]
                    rec["seed_ok"] = bool(r["service_ok"] or r["n_fcs_ok"])
                    if "n_codewords" in r:
                        rec["ldpc_codewords_ok"] = f"{r['n_codewords_ok']}/{r['n_codewords']}"
                    macs = [s["mac"] for s in r["subframes"] if s.get("mac")]
                    if macs:
                        rec["mac"], rec["macs"] = macs[0], macs
                except Exception as exc:                # noqa: BLE001 — record, don't fail the run
                    rec["payload_error"] = str(exc)[:120]
        elif p.fmt == "VHT" and p.sig_crc_ok:
            # The payload is out of reach, but VHT-SIG-A (on the primary 20 MHz) still says who sent
            # what to whom: Group ID gives the direction, NSS/LDPC the device class.
            e = p.extra
            gid = e.get("group_id")
            rec["nss"], rec["ldpc"] = e.get("nsts_su"), bool(e.get("ldpc"))
            rec["direction"] = VHT_DIRECTION.get(gid, "multi-user MIMO" if gid is not None else None)
            if e.get("partial_aid") is not None:
                rec["partial_aid"] = e["partial_aid"]
            rec["payload_note"] = _vht_note(p.bandwidth_mhz, e.get("nsts_su"), e.get("ldpc"), cap.fs)
        elif p.fmt == "HE" and p.sig_crc_ok:
            # HE-SIG-A (CRC-4): UL/DL bit, BSS colour (ties a PPDU to its AP), MCS / streams / coding / TXOP
            e = p.extra
            rec.update(he_format=e["he_format"], bss_color=e["bss_color"], txop_us=e.get("txop_us"),
                       direction=HE_DIRECTION[bool(e["uplink"])])
            if e.get("nsts") is not None:
                rec.update(nss=e["nsts"], ldpc=bool(e.get("ldpc")), dcm=bool(e.get("dcm")))
            rec["payload_note"] = _he_note(p.bandwidth_mhz, e.get("nsts"), cap.fs)
        out.append(rec)
    return out


def _scan(ctx: DecodeContext):
    if "wifi_ofdm" in ctx.cache:
        return ctx.cache["wifi_ofdm"]
    from . import ofdm_rx as W
    cap = ctx.cap
    lo, hi = cap.fc - cap.fs / 2, cap.fc + cap.fs / 2
    det = W.WifiDetector(cap.fs, cap.fc)
    pk = det.run(cap.iq) if WIFI_BANDS & set(band_keys(lo, hi)) else []
    n_rej = len(det.rejects)
    # 802.11 off the standard channel grid / outside the Wi-Fi bands (fixed-wireless links, telemetry radios):
    # blind 1 MHz centres only where a 3.2 us CP correlation remains after the known PPDUs are blanked
    extra = W.offgrid_plan(cap.iq, cap.fs, cap.fc, pk, X=ctx.ch.X)
    if extra:
        d2 = W.WifiDetector(cap.fs, cap.fc, channels=extra)
        pk = W.merge_packets(pk, d2.run(cap.iq))
        n_rej += len(d2.rejects)
    out = _records(pk, cap)
    _attribute_vht(out)
    ctx.cache["wifi_ofdm"] = (out, n_rej)
    return ctx.cache["wifi_ofdm"]


def _scan_scaled(ctx: DecodeContext):
    """[(link, records)] for 802.11 links at a non-standard clock / frequency / orientation."""
    if "wifi_ofdm_scaled" in ctx.cache:
        return ctx.cache["wifi_ofdm_scaled"]
    from . import ofdm_scaled as S
    cap = ctx.cap
    out = [(ln, _records(ln.packets, cap, ln.clock_hz)) for ln in S.scan(cap.iq, cap.fs, cap.fc)]
    for ln, recs in out:
        ln.packets = []                                 # the records hold what is needed; free the IQ state
    ctx.cache["wifi_ofdm_scaled"] = out
    return out


VHT_DIRECTION = {0: "uplink (station to AP)", 63: "downlink (AP to station)"}


def _block_s(ctx):
    return len(ctx.cap.iq) / ctx.cap.fs


def _vht_note(bw_mhz, nss, ldpc, fs):
    """Why this VHT payload is not decoded -- the first reason that applies, said plainly."""
    if bw_mhz and bw_mhz * 1e6 > fs:
        return f"{bw_mhz} MHz payload is wider than the {fs / 1e6:g} MHz capture"
    if nss and nss > 1:
        return f"{nss} spatial streams (needs as many antennas)"
    if ldpc:
        return "LDPC-coded VHT payload (only 20 MHz HT LDPC implemented)"
    return f"{bw_mhz} MHz VHT payload decoding not implemented"


HE_DIRECTION = {True: "uplink (station to AP)", False: "downlink (AP to station)"}


def _he_note(bw_mhz, nss, fs):
    """Why this HE payload is not decoded."""
    if bw_mhz and bw_mhz * 1e6 > fs:
        return f"{bw_mhz} MHz payload is wider than the {fs / 1e6:g} MHz capture"
    if nss and nss > 1:
        return f"{nss} spatial streams (needs as many antennas)"
    return "HE (802.11ax) payload decoding not implemented"


def _he_fields(mine, all_recs):
    """HE summary: airtime by direction, BSS colours, and colour -> BSSID from the beacons' HE Operation
    element (only when exactly one BSSID in the capture announces that colour)."""
    he = [p for p in mine if p["fmt"] == "HE" and p.get("bss_color") is not None]
    if not he:
        return {}
    air = Counter()
    for p in he:
        air[p["direction"]] += p.get("dur_us") or 0
    col = {}
    for r in all_recs:
        for m in r.get("macs", []):
            c = (m.get("he_operation") or {}).get("bss_color")
            if c and m.get("bssid"):
                col.setdefault(c, set()).add(m["bssid"])
    f = {"he_airtime_ms_by_direction": {k: round(v / 1e3, 3) for k, v in air.most_common()},
         "he_formats": dict(Counter(p["he_format"] for p in he)),
         "he_bss_colors": {str(k): v for k, v in sorted(Counter(p["bss_color"] for p in he).items())}}
    hits = {str(c): sorted(col[c])[0] for c in sorted({p["bss_color"] for p in he}) if len(col.get(c, ())) == 1}
    if hits:
        f["he_bss_from_color"] = hits
    return f


def _attribute_vht(recs):
    """Tie uplink VHT PPDUs to an access point through the Partial AID (802.11-2020 26.17).

    A station sending to its AP puts nine bits of the BSSID in VHT-SIG-A. Many clients set it to 0
    instead (seen on every uplink in the 2026-10-05 live 5 GHz captures), and 0 is also what an
    unknown value looks like, so 0 is never matched. A match is reported only when the value is
    non-zero AND exactly one BSSID decoded in this capture produces it.
    """
    from . import payload as WP

    bssids = sorted({m.get("bssid") for r in recs for m in r.get("macs", [])
                     if m.get("bssid") and not int(m["bssid"][:2], 16) & 1})
    paid = {}
    for b in bssids:
        paid.setdefault(WP.vht_uplink_paid(b), []).append(b)
    for r in recs:
        if (r.get("direction") or "").startswith("uplink") and r.get("partial_aid"):
            hit = paid.get(r["partial_aid"], [])
            if len(hit) == 1:
                r["bssid_from_partial_aid"] = hit[0]


STEPS_COMMON = [
    "non-HT PPDUs: OFDM demod, deinterleave, Viterbi, descramble, CRC-32 FCS (receivers.wifi.payload)",
    "HT (802.11n) MCS 0-7, 20 MHz: channel per space-time stream from the HT-LTFs (smoothed when the HT-SIG "
    "allows), Alamouti combining for STBC, pilots p_3.. with the rotating HT pattern, per-symbol phase and "
    "timing-drift tracking; BCC: N_COL=13 deinterleave + Viterbi; LDPC: 802.11 codeword layout "
    "(shortening/puncturing/repetition), log-MAP LLRs, min-sum decoding, data-aided channel re-estimation "
    "from codewords that pass their parity checks; descramble; an A-MPDU is split on its delimiters and "
    "every MPDU checked on its own FCS",
    "MAC dissection (mac80211.py) incl. Remote ID (ASTM F3411 / DJI) when a drone broadcasts it",
    "transmitter attribution (links.py): MAC header, control responses (ACK/BA/CTS RA), sequential "
    "scrambler-seed counters read from the SERVICE field; received power as a weak hint"]


def _field_frames(recs):
    """PPDU records -> field.py frames: airtime, PSDU bytes when known, transmitter, MPDUs, control / video flags."""
    out = []
    for p in recs:
        ms = p.get("macs") or ([p["mac"]] if p.get("mac") else [])
        b = p["length"] if p["fmt"] == "NON_HT" else p.get("psdu_bytes") if p["fmt"] == "HT" else None
        ctrl = any(m.get("type") == "ctrl" for m in ms) or (not ms and b is not None and b <= CTRL_MAX_BYTES)
        out.append({"t0": p["t0"], "t1": p["t1"], "dur_us": p.get("dur_us") or 0.0, "bytes": b, "ta": p.get("ta"),
                    "macs": ms, "ctrl": ctrl,
                    "video": any((m.get("ipv4") or {}).get("video") for m in ms)})
    return out


def field_summary(frames, block_s):
    """BSS link profiles and drone-link verdicts of one signal (field.py). Shared by the OFDM and DSSS plugins."""
    from . import field as FD
    out = {}
    lk = FD.links(frames, block_s)
    if lk:
        out["wifi_links"] = lk[:10]
        dl = [{"bssid": r["bssid"], "ssid": r.get("ssid"), "clients": r["clients"], **r["drone_link"]}
              for r in lk if r.get("drone_link")]
        if dl:
            out["drone_links"] = dl
    return out


def drone_label(fields):
    """Label suffix and evidence for a signal carrying a drone link."""
    dl = fields.get("drone_links")
    if not dl:
        return "", ""
    d = max(dl, key=lambda x: x["confidence"] or 0)
    what = d["platform"] or "by behaviour"
    return (f", drone link-like ({what})",
            f"; drone link-like on BSS {d['bssid']} ({what}, {d['basis']}, confidence {d['confidence']})")


CTRL_MAX_BYTES = 40


def _fields(mine, sifs_us, block_s=None):
    """Shared summary fields: counts, MACs, attribution."""
    from . import links
    fmts = Counter(p["fmt"] for p in mine)
    fcs = [p for p in mine if p.get("fcs_ok")]
    macs = [m for p in fcs for m in p.get("macs", [p["mac"]] if "mac" in p else [])]
    tx = links.attribute(mine, sifs_us=sifs_us)
    f = {"n_ppdus": len(mine), "formats": dict(fmts), "channels": dict(Counter(p["channel"] for p in mine)),
         "bandwidths_mhz": dict(Counter(p["bw_mhz"] for p in mine)), "n_fcs_ok": len(fcs), "n_mpdus": len(macs),
         "bssids": sorted({m.get("bssid") for m in macs if m.get("bssid")}),
         "ssids": sorted({m["ssid"] for m in macs if m.get("ssid")}),
         "transmitters": tx,
         "n_attributed": sum(1 for p in mine if p.get("ta")),
         "ppdus": [{k: v for k, v in p.items() if k not in ("f0", "f1", "t0", "t1")} for p in mine[:200]]}
    n_ht = sum(1 for p in mine if p["fmt"] == "HT")
    if n_ht:
        f["n_ht_ppdus"] = n_ht
        f["n_ht_decoded"] = sum(1 for p in mine if p["fmt"] == "HT" and p.get("fcs_ok"))
    rids = [m["remote_id"] for m in macs if m.get("remote_id")]
    if rids:                                         # research, PLAN M2: manufacturer from the Remote ID serial
        from spectrum_intel.core.drone_models import remote_id_model
        f["model"] = remote_id_model(rids[-1])
    notes = Counter(p["payload_note"] for p in mine if p.get("payload_note"))
    if notes:
        f["payload_not_decodable"] = dict(notes)
    if block_s:                                      # research, PLAN W1
        f.update(field_summary(_field_frames(mine), block_s))
    return f, fcs


class WifiOfdmDecoder(Decoder):
    name = "wifi_ofdm"
    protocol = "IEEE 802.11 OFDM"
    capture_level = True

    def can_decode(self, sig, ctx):
        return bool(WIFI_BANDS & set(band_keys(sig.f0, sig.f1))) and sig.bw >= 2e6

    def capture_signals(self, ctx):
        pk, _ = _scan(ctx)
        groups = {}
        for p in pk:
            groups.setdefault((p["channel"], p["bw_mhz"]), []).append(p)
        out = []
        for (chn, bw), g in groups.items():
            res = self.summarise(g, len(pk), 0, pk, _block_s(ctx))
            f0, f1 = res["span"]
            snr = float(np.median([p["snr_db"] for p in g]))
            out.append(receiver_signal([(p["t0"], p["t1"]) for p in g], f0, f1, res, snr, self.name))
        return out

    def decode(self, sig, ctx):
        pk, n_rej = _scan(ctx)
        mine = [p for p in pk if overlaps(sig, p["t0"], p["t1"], p["f0"], p["f1"])]
        return self.summarise(mine, len(pk), n_rej, pk, _block_s(ctx))

    def summarise(self, mine, n_all, n_rej, all_recs=(), block_s=None):
        steps = ["capture-wide 802.11 legacy-preamble scan (receivers.wifi.ofdm_rx.WifiDetector.run): every 20 MHz "
                 "channel resampled to 20 MSps, L-STF autocorrelation plateau, L-LTF channel estimate, "
                 "L-SIG Viterbi + parity + bit-agreement check, HT-SIG/VHT-SIG-A CRC-8, HE-SIG-A CRC-4 (BSS colour, UL/DL, MCS, streams, TXOP)",
                 "off the standard grid / outside the Wi-Fi bands: 20 MHz slices every 5 MHz with the decoded PPDUs "
                 "blanked; >= 3 adjacent slices with a 3.2 us CP correlation -> blind 1 MHz centres around them "
                 "(ofdm_rx.offgrid_plan)",
                 f"{n_all} PPDUs in capture, {n_rej} preambles rejected; {len(mine)} attributed to this signal"
                 ] + STEPS_COMMON
        if not mine:
            return result("not_possible", "receivers.wifi ofdm_rx L-STF/L-LTF/L-SIG receiver", steps,
                          reason="no 802.11 OFDM preamble found in this signal's time-frequency span")
        fields, fcs = _fields(mine, 16.0, block_s)
        fmts = fields["formats"]
        vht = [p for p in mine if p.get("direction") and p["fmt"] == "VHT"]
        if vht:
            air = Counter()
            for p in vht:
                air[p["direction"]] += p.get("dur_us") or 0
            fields["vht_airtime_ms_by_direction"] = {k: round(v / 1e3, 3) for k, v in air.most_common()}
            fields["vht_streams"] = dict(Counter(p.get("nss") for p in vht))
            hits = Counter(p["bssid_from_partial_aid"] for p in vht if p.get("bssid_from_partial_aid"))
            if hits:
                fields["vht_uplink_bss_from_partial_aid"] = dict(hits)
        fields.update(_he_fields(mine, all_recs or mine))
        status = "success" if fcs else "partial"
        reason = ("" if fcs else "PHY headers (L-SIG/HT-SIG) decoded; payloads need multiple antennas "
                  "(2+ spatial streams) or a wider channel, or they failed FCS")
        ch = Counter(p["channel"] for p in mine).most_common(1)[0][0]
        suf, dev = drone_label(fields)
        lab = f"Wi-Fi 802.11 OFDM ({'/'.join(sorted(fmts))}) {ch}{suf}"
        ev = (f"{len(mine)} 802.11 PPDUs decoded (L-SIG parity/agreement" +
              (", HT/VHT/HE-SIG CRC" if any(p['sig_crc_ok'] for p in mine) else "") + ")" +
              (f"; {len(fcs)} MPDUs pass CRC-32 FCS" if fcs else "") + dev)
        f_lo = min(p["f0"] for p in mine)
        f_hi = max(p["f1"] for p in mine)
        return result(status, "receivers.wifi ofdm_rx + payload (research receiver)", steps, fields,
                      reason, label=lab, confidence=0.97 if fcs else 0.9,
                      modulation="OFDM 64-pt, 312.5 kHz SCS (BPSK..1024-QAM)", evidence=ev, span=(f_lo, f_hi))


class WifiOfdmScaledDecoder(Decoder):
    """802.11 at a non-standard clock, frequency or spectrum orientation (any band)."""
    name = "wifi_ofdm_scaled"
    protocol = "IEEE 802.11 OFDM, non-standard numerology"
    capture_level = True

    def capture_signals(self, ctx):
        out = []
        for ln, recs in _scan_scaled(ctx):
            if not recs:
                continue
            res = self.summarise(ln, recs, _block_s(ctx))
            f0, f1 = res["span"]
            snr = float(np.median([p["snr_db"] for p in recs]))
            out.append(receiver_signal([(p["t0"], p["t1"]) for p in recs], f0, f1, res, snr, self.name))
        return out

    def summarise(self, ln, recs, block_s=None):
        nm = ln.numerology
        clk = ln.clock_hz / 1e6
        steps = ["bursty wideband bands with repeatable edges from a ~60 kHz spectrogram (ofdm_scaled.find_bands)",
                 f"numerology of the band: CP autocorrelation peak at T_u = {nm['tu_s'] * 1e6:.3f} us "
                 f"({nm['cp_corr']:.2f}, {nm['cp_corr_over_median']}x median), symbol T_s = {nm['ts_s'] * 1e6:.3f} us, "
                 f"guard {nm['guard']} of T_u (802.11: 0.25) -> 64-subcarrier clock {nm['clock_measured_hz'] / 1e6:.4f} "
                 f"MHz, taken as {clk:g} MHz; occupied {nm['occupied_hz'] / 1e6:.2f} MHz = "
                 f"{nm['occupied_over_clock']} x clock",
                 f"capture time-scaled by 20/{clk:g} into a standard 20 MHz channel"
                 + (" after mirroring the spectrum (complex conjugate): only the mirrored orientation gives "
                    "parity/CRC-valid SIG fields" if ln.inverted else "")
                 + "; ofdm_rx: L-STF plateau, L-LTF, L-SIG parity + bit agreement, HT-SIG CRC",
                 f"{len(recs)} PPDUs, {ln.n_rejects} preambles rejected"] + STEPS_COMMON
        fields, fcs = _fields(recs, 16.0 * 20 / clk, block_s)
        fmts = fields["formats"]
        fields.update(center_mhz=round(ln.center_hz / 1e6, 4), clock_mhz=clk, channel_width_mhz=clk,
                      subcarrier_spacing_khz=round(312.5 * clk / 20, 3), spectrum_inverted=ln.inverted,
                      numerology={k: v for k, v in nm.items() if k not in ("reject",)},
                      standard_equivalent=("802.11 quarter-clocked (5 MHz) channel" if clk == 5 else
                                           "802.11 half-clocked (10 MHz) channel" if clk == 10 else
                                           "802.11 at 20 MHz, off the channel grid" if clk == 20 else
                                           f"802.11 run at a {clk:g} MHz clock (non-standard)"))
        vendors = sorted({t["vendor"] for t in fields["transmitters"]
                          if t.get("vendor") and "unknown" not in t["vendor"].lower()
                          and t["evidence"].get("mac header")})
        phy = "802.11n" if "HT" in fmts else "802.11a/g"
        lab = (f"Wi-Fi {phy}-type link, {clk:g} MHz channel at {ln.center_hz / 1e6:.2f} MHz"
               + (", spectrum-inverted" if ln.inverted else "") + (f" ({', '.join(vendors)})" if vendors else ""))
        suf, dev = drone_label(fields)
        lab += suf
        ev = (f"{len(recs)} PPDUs with parity/CRC-checked SIG fields at a {clk:g} MHz OFDM clock"
              + (f"; {len(fcs)} PPDUs carry CRC-32-valid MPDUs" if fcs else "")
              + (f"; transmitters {', '.join(t['ta'] for t in fields['transmitters'])}" if fields["transmitters"] else "")
              + dev)
        status = "success" if fcs else "partial"
        reason = "" if fcs else "PHY headers decoded; no MPDU passed its FCS"
        f_lo = min(p["f0"] for p in recs)
        f_hi = max(p["f1"] for p in recs)
        return result(status, "receivers.wifi ofdm_scaled + ofdm_rx + payload (time-scaled 802.11)", steps, fields,
                      reason, label=lab, confidence=0.97 if fcs else 0.9,
                      modulation=f"OFDM 64-pt, {312.5 * clk / 20:.2f} kHz SCS (802.11 at a {clk:g} MHz clock)",
                      evidence=ev, span=(f_lo, f_hi))


PLUGINS = [WifiOfdmDecoder(), WifiOfdmScaledDecoder()]
