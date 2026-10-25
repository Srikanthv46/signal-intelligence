"""Pipeline plugin: LoRa/ELRS, FSK RC links, BLE adverts and 2 Mb/s ESB-like links (wraps .detector)."""
from __future__ import annotations


from spectrum_intel.core.bands_india import BLE_ADV
from spectrum_intel.decode.emission_adapter import EmissionReceiver

RC_BANDS = {"ism_24", "ism_433", "ism_865", "b8_gap"}
# sub-GHz LoRa links outside the delicensed bands: EU868, FCC915 drone links (ELRS / Crossfire 902-928 MHz) and
# off-band telemetry radios (live: SX126x SF5 at 890.03 MHz). A capture that touches no RC band but this range
# gets the LoRa scan only (preamble + SFD verified; FSK would mistake GSM/LTE uplink bursts for RC links).
LORA_HZ = (860e6, 930e6)


def _lora_overlap(lo, hi):
    return lo < LORA_HZ[1] and hi > LORA_HZ[0]
ESB_LABEL = "2.4 GHz GFSK link (ESB-like, 2 Mb/s)"
ADV_TOL_HZ = 0.35e6          # BLE advertising channel centre +- (receiver ppm + carrier offset), channels 2 MHz apart


class RcLinkDecoder(EmissionReceiver):
    name = "rc_link"
    protocol = "LoRa/FSK/BLE/ESB"
    key = "rc_link"
    bands = RC_BANDS
    method = "receivers.rc_link (research receiver)"
    steps = ["LoRa: blind dechirp scan per sub-band (critically sampled, K-symbol coherent FFT), 2.25-downchirp "
             "SFD verification, symbol walk -> SF/BW/preamble/payload length, ELRS channel-grid match",
             "GFSK/FSK: matched-area energy at 3 widths, narrowband flank test, constant-envelope test, "
             "inst.-frequency clusters (deviation), symbol rate from transition timing",
             "1 Mb/s h~0.5 bursts: BLE preamble + access address 0x8E89BED6, de-whitening, CRC-24",
             "2.4 GHz 0.5-3 MHz bursts: 2 Mb/s GFSK slicing, >= 16-bit 0101 preamble -> ESB-like packet; packets "
             "sharing a sync word: CRC-16/CCITT end and init found jointly (>= 3 distinct packets on one CRC)"]

    def _rc_capture(self, ctx):
        return super().in_range(ctx.cap.fc - ctx.cap.fs / 2, ctx.cap.fc + ctx.cap.fs / 2)

    def in_range(self, lo, hi):
        return super().in_range(lo, hi) or _lora_overlap(lo, hi)

    def can_decode(self, sig, ctx):
        if self._rc_capture(ctx):
            return super().can_decode(sig, ctx)
        lo, hi = ctx.cap.fc - ctx.cap.fs / 2, ctx.cap.fc + ctx.cap.fs / 2
        return _lora_overlap(lo, hi) and _lora_overlap(sig.f0, sig.f1)

    def run_detector(self, ctx):
        from .detector import detect
        fc = ctx.cap.fc
        lora_only = None if self._rc_capture(ctx) else (LORA_HZ[0] - fc, LORA_HZ[1] - fc)
        return detect(ctx.cap.iq, ctx.cap.fs, fc, {"lora_only": lora_only} if lora_only else {})

    def group_key(self, e):
        if e.family == "esb_link":                   # ESB channels sit on the 1 MHz raster (2400 + k MHz)
            return e.label, round(e.center_mhz)
        return super().group_key(e)

    def pretty(self, top, em):
        if any(e.family == "esb_link" for e in em):
            return ESB_LABEL
        # by family, not substring: "lora_preamble" contains "ble"
        if any(e.family == "bluetooth" for e in em) or any((e.decoded or {}).get("aa") for e in em):
            d = next((e.decoded for e in em if (e.decoded or {}).get("crc_ok")), None)
            if d:
                lr = " Long Range" if str(d.get("phy", "")).startswith("LE Coded") else ""
                if any(((e.decoded or {}).get("ad") or {}).get("remote_id") for e in em):
                    return f"Drone Remote ID (Bluetooth LE{lr}, {d.get('pdu_name')})"
                trk = next((dv for e in em if (e.decoded or {}).get("crc_ok")
                            for dv in [((e.decoded or {}).get("ad") or {}).get("device")] if dv and dv["tracker"]), None)
                if trk:                                  # item tracker (research, signal gaps S7)
                    return f"Bluetooth LE tracker ({trk['label']})"
                return f"Bluetooth LE{lr} advertising ({d.get('pdu_name')})"
            # advertising uses only channels 37/38/39; a 1 Mb/s GFSK hop elsewhere is a BLE data channel or a
            # BR/EDR hop (the measured modulation index does not separate h 0.5 from 0.32 reliably)
            br = next((e.decoded for e in em if (e.decoded or {}).get("lap")), None)
            if br:                                       # BCH-verified BR/EDR access code (research, signal gaps S7)
                return f"Bluetooth BR/EDR (LAP {br['lap']}{', inquiry' if br['inquiry'] else ''})"
            if any(abs(e.center_mhz * 1e6 - f) <= ADV_TOL_HZ for e in em for f in BLE_ADV.values()):
                return "Bluetooth LE advertising-like"
            return "Bluetooth (BLE/BR)-like"
        if any(e.label == "sik_telemetry" for e in em):
            return "SiK telemetry radio (MAVLink)"
        fr = [e.decoded for e in em if e.label.startswith("frsky_") and e.decoded]
        if fr:
            d = next((d for d in fr if d.get("crc_ok")), fr[0])
            return d["protocol"]
        zb = [e.decoded for e in em if e.label.startswith("ieee802154") and e.decoded]
        if zb:
            ok = [d for d in zb if d.get("crc_ok")]
            v = next((d["vendor"] for d in ok if d.get("vendor")), None)
            return (f"IEEE 802.15.4 frames ({v})" if v else "IEEE 802.15.4 frames (XBee / Zigbee-like)") if ok \
                else "IEEE 802.15.4-like (O-QPSK 2.4 GHz)"
        ty = [e.decoded for e in em if e.label == "toy_drone" and e.decoded]
        if ty:
            return f"Toy drone RC link: {ty[0]['protocol']}"
        xn = [e.decoded for e in em if e.label in ("bayang", "xn297_like") and e.decoded]
        if xn:
            return next((d for d in xn if d.get("checksum_ok")), xn[0])["protocol"]
        af = [e.decoded for e in em if e.label.startswith(("afhds2a", "afhds1")) and e.decoded]
        if af:
            return next((d for d in af if d.get("crc_ok")), af[0])["protocol"]
        if any(e.label == "elrs_gfsk" for e in em):
            return "ExpressLRS GFSK (LR1121)"
        rl = [e.decoded for e in em if e.label == "radiolink" and e.decoded]
        if rl:
            return next((d for d in rl if d.get("packet_name") == "data"), rl[0])["protocol"]
        hb = [e.decoded for e in em if e.label == "hubsan" and e.decoded]
        if hb:
            return "Hubsan (A7105)"
        cc = [e.decoded for e in em if e.label in ("hott", "hitec", "sfhss") and e.decoded]
        if cc:
            return cc[0]["protocol"].replace(" receiver telemetry", "")
        og = [e.decoded for e in em if e.label == "ogn_like" and e.decoded]
        if og:
            return "OGN tracker-like (OGNTP)" + (", UAV" if any(d.get("uav") for d in og) else "")
        ad = [e.decoded for e in em if e.label == "adsl" and e.decoded]
        if ad:
            return "ADS-L iConspicuity (SRD-860)" + (", UAS" if any(d.get("uas") for d in ad) else "")
        fa = [e.decoded for e in em if (e.decoded or {}).get("protocol") == "FANET" and e.decoded.get("crc_ok")]
        if fa:
            uav = any((d.get("fields") or {}).get("uav") for d in fa)
            return "FANET aircraft tracking (LoRa)" + (", UAV" if uav else "")
        if any(e.label == "flarm_like" for e in em):
            return "FLARM-like (868 MHz collision-avoidance beacon)"
        if any(e.label == "crossfire_sync" for e in em):
            return "TBS Crossfire-like (150 Hz FSK, sync 0x3E6B3E)"
        if "mlrs" in top:
            return top.replace("_", " ").replace("mlrs", "mLRS").replace("lora", "LoRa")
        if "lora" in top or "elrs" in top:
            return top.replace("_", " ").replace("elrs", "ExpressLRS").replace("lora", "LoRa")
        return top.replace("_", " ")

    def summarise(self, em, sig):
        out = super().summarise(em, sig)
        elrs = [e for e in em if "elrs" in e.label and e.features.get("modulation") == "lora"
                and e.center_mhz >= 2400.0 and e.features.get("iq_inverted") is not None]
        if elrs:
            # ExpressLRS 2.4 GHz (SX1280): IQ inverted <=> UID[5] bit 0 = 1, or binding mode (always inverted); the
            # SX127x sub-GHz driver never inverts (research: research/drone_transmitters/07_elrs_details.md)
            inv = {bool(e.features["iq_inverted"]) for e in elrs}
            if len(inv) == 1:
                v = inv.pop()
                out["fields"]["elrs_identity"] = {
                    "uid5_bit0": int(v), "iq_inverted": v, "packets": len(elrs),
                    "note": "binding phrase UID[5] bit 0 from the LoRa IQ polarity (1 also while binding)"}
        syncs = [e.decoded for e in em if (e.decoded or {}).get("protocol") == "ExpressLRS"
                 and e.decoded.get("packet_type") == "SYNC" and e.decoded.get("crc_ok")]
        if syncs:
            # ExpressLRS SYNC packets verified by their own UID bytes (research, PLAN R2c)
            d = syncs[-1]
            f = {k: v for k, v in d["sync"].items()}
            f.update(packets=len(syncs), ota_version=d["ota_version"], binding=d["binding"])
            if d.get("hop_recovery"):                    # research, PLAN R2d
                f["hop_recovery"] = d["hop_recovery"]
            out["fields"]["elrs_sync"] = f
        og = [e.decoded for e in em if e.label == "ogn_like" and e.decoded]
        if og:
            # OGN tracker (research, PLAN G10): addresses, aircraft type (13 = UAV), last position; no LDPC check
            pos = [d for d in og if d.get("lat") is not None]
            out["fields"]["ogn"] = {
                "packets": len(og), "addresses": sorted({f"{d['address_type']}:{d['address']}" for d in og}),
                "aircraft_types": sorted({d["aircraft_type"] for d in og if d.get("aircraft_type")}),
                "uav": any(d.get("uav") for d in og), "encrypted": any(d["encrypted"] for d in og),
                "position": {k: pos[-1][k] for k in ("lat", "lon", "alt_m")} if pos else None,
                "fec": "LDPC not checked"}
        ad = [e.decoded for e in em if e.label == "adsl" and e.decoded]
        if ad:
            # ADS-L SRD-860 (research, PLAN G10): addresses, payload types, category (11-13 UAS), last position
            pos = [d for d in ad if d.get("lat") is not None]
            out["fields"]["adsl"] = {
                "packets": len(ad), "addresses": sorted({f"{d['address_table']}:{d['address']}" for d in ad}),
                "payload_types": sorted({d["payload_type"] for d in ad}),
                "aircraft_categories": sorted({d["aircraft_category"] for d in ad if d.get("aircraft_category")}),
                "uas": any(d.get("uas") for d in ad),
                "position": {k: pos[-1][k] for k in ("lat", "lon", "alt_m") if pos[-1].get(k) is not None} if pos else None}
        fa = [e.decoded["fields"] for e in em if (e.decoded or {}).get("protocol") == "FANET"
              and e.decoded.get("crc_ok") and e.decoded.get("fields")]
        if fa:
            # FANET (research, PLAN G10): device IDs, aircraft type (7 = UAV), last position, names
            pos = [d for d in fa if d.get("lat") is not None]
            out["fields"]["fanet"] = {
                "packets": len(fa), "device_ids": sorted({d["device_id"] for d in fa}),
                "manufacturers": sorted({d["manufacturer"] for d in fa}),
                "message_types": sorted({d["message_type"] for d in fa}),
                "aircraft_types": sorted({d["aircraft_type"] for d in fa if d.get("aircraft_type")}),
                "uav": any(d.get("uav") for d in fa),
                "position": {k: pos[-1][k] for k in ("lat", "lon", "alt_m") if pos[-1].get(k) is not None} if pos else None,
                "names": sorted({d["name"] for d in fa if d.get("name")})}
        fsk = [e for e in em if e.family in ("rc_fhss", "narrowband", "unknown")
               and e.features.get("symbol_rate") and e.features.get("modulation") == "gfsk"]
        if fsk:
            # research catalog of RC / telemetry protocols (rate, deviation, band, packet period)
            from .protocols import candidates
            import numpy as np
            R = float(np.median([e.features["symbol_rate"] for e in fsk]))
            dv = [e.features["dev_hz"] for e in fsk if e.features.get("dev_hz")]
            st = sorted(e.start_us for e in fsk)
            gaps = [b - a for a, b in zip(st[:-1], st[1:]) if b - a > 300.0]
            iv = min(gaps) / 1e3 if len(st) >= 3 and gaps else None
            pc = candidates(R, float(np.median(dv)) if dv else None, float(np.median([e.center_mhz for e in fsk])), iv)
            if pc:
                out["fields"]["protocol_candidates"] = pc
        rids = [((e.decoded or {}).get("ad") or {}).get("remote_id") for e in em if (e.decoded or {}).get("crc_ok")]
        rids = [r for r in rids if r]
        if rids:
            # Remote ID (F3411 over BLE; research capability): later messages complete / update earlier ones
            merged: dict = {}
            for r in rids:
                for k, v in r.items():
                    if k == "messages":
                        merged.setdefault("messages", [])
                        merged["messages"] += [m for m in v if m not in merged["messages"]]
                    elif v is not None:
                        merged[k] = v
            f = out["fields"]
            f["remote_id"] = merged
            f["remote_id_adverts"] = len(rids)
            if merged.get("uas_id"):
                f["uas_id"] = merged["uas_id"]
            from spectrum_intel.core.drone_models import remote_id_model
            f["model"] = remote_id_model(merged)                       # research, PLAN M2
            out["evidence"] = (out.get("evidence") or "") + (f"; Remote ID (ASTM F3411) in {len(rids)} CRC-valid advert(s): "
                                                             f"{', '.join(merged.get('messages', []))}")
        sik = [e.decoded for e in em if e.label == "sik_telemetry" and (e.decoded or {}).get("crc_ok")]
        if sik:
            # SiK telemetry (research): radio facts + the vehicle state from the CRC-valid MAVLink frames
            from spectrum_intel.core import mavlink as MAV
            frames = [f for d in sik for f in d.get("mavlink_frames", [])]
            out["fields"]["sik"] = {"netid": sik[0]["netid"], "air_rate_bps": sik[0]["air_rate_bps"],
                                    "whitening": sik[0]["whitening"], "hw_crc_variant": sik[0]["hw_crc_variant"],
                                    "packets": len(sik), "mavlink_frames_crc_ok": sum(d["mavlink_crc_ok"] for d in sik)}
            mv = MAV.summary(frames)
            if mv:
                out["fields"]["mavlink"] = mv
                out["evidence"] = (out.get("evidence") or "") + (
                    f"; SiK NETID {sik[0]['netid']}, MAVLink {', '.join(mv.get('messages', []))}")
        ads = [(e.decoded or {}).get("ad") for e in em if (e.decoded or {}).get("ad")]
        if ads:
            keep = {}
            for a in ads:
                for k in ("local_name", "company_id", "company", "service_uuids16", "tx_power_dbm", "appearance"):
                    if a.get(k) is not None:
                        keep[k] = a[k]
            if keep:
                out["fields"]["ble_ad"] = keep
            # trackers / notable devices (research, signal gaps S7): distinct recognitions of CRC-valid adverts
            devs, seen = [], set()
            for e in em:
                dv = ((e.decoded or {}).get("ad") or {}).get("device") if (e.decoded or {}).get("crc_ok") else None
                if dv:
                    k = (dv["label"], e.decoded.get("adv_addr"))
                    if k not in seen:
                        seen.add(k)
                        devs.append({**dv, "adv_addr": e.decoded.get("adv_addr")})
            if devs:
                out["fields"]["ble_devices"] = devs[:32]
                trk = sorted({d["label"] for d in devs if d["tracker"]})
                if trk:
                    out["evidence"] = (out.get("evidence") or "") + f"; item tracker(s) in the adverts: {', '.join(trk)}"
        br = [e.decoded for e in em if e.label == "bluetooth_bredr" and e.decoded]
        if br:
            # Bluetooth BR/EDR (research, signal gaps S7): LAPs of the BCH-verified access codes; UAP / NAP not on air
            out["fields"]["bredr"] = {"laps": sorted({d["lap"] for d in br}), "packets": len(br),
                                      "inquiry": any(d["inquiry"] for d in br),
                                      "max_access_code_errors": max(d["access_code_errors"] for d in br)}
            out["evidence"] = (out.get("evidence") or "") + (
                f"; {len(br)} BR/EDR access code(s) within {out['fields']['bredr']['max_access_code_errors']} bit(s) "
                f"of the (64,30) BCH codeword of LAP {', '.join(out['fields']['bredr']['laps'])}")
        cf = [e.decoded for e in em if e.label == "crossfire_sync" and e.decoded]
        if cf:
            # TBS Crossfire 150 Hz (research, PLAN R9): sync-word packets, grid channels and links seen
            grid = [d["grid"] for d in cf if d.get("grid")]
            out["fields"]["crossfire"] = {
                "packets": len(cf), "sync": cf[0]["sync"],
                "band": grid[0]["band"] if grid else None,
                "channels": sorted({(g["link"], g["channel"]) for g in grid}),
                "note": "packets read but not CRC-verified (per-hop CRC init unknown)"}
        fr = [e.decoded for e in em if e.label.startswith("frsky_") and e.decoded]
        if fr:
            # FrSky D8 / D16 (research, PLAN R6): variant from the CRC structure + rate, TX ID from bytes 1-2
            ok = [d for d in fr if d.get("crc_ok")]
            out["fields"]["frsky"] = {
                "packets": len(fr), "crc_ok": len(ok), "variants": sorted({d["variant"] for d in (ok or fr)}),
                "tx_ids": sorted({d["tx_id"] for d in ok if d.get("tx_id")}), "sync": fr[0]["sync"]}
        af = [e.decoded for e in em if e.label.startswith("afhds2a") and e.decoded]
        if af:
            # FlySky AFHDS 2A (research, PLAN R7): ID-code packets; CRC-verified ones give the IDs and packet types
            ok = [d for d in af if d.get("crc_ok")]
            out["fields"]["afhds2a"] = {
                "packets": len(af), "crc_ok": len(ok), "tx_ids": sorted({d["tx_id"] for d in ok}),
                "rx_ids": sorted({d["rx_id"] for d in ok}), "packet_types": sorted({d["packet_type"] for d in ok}),
                "binding": any(d["binding"] for d in ok),
                "convention": ok[0]["convention"] if ok else None}
        a1 = [e.decoded for e in em if e.label.startswith("afhds1") and e.decoded]
        if a1:
            # FlySky AFHDS 1 (research, PLAN G5): TX ID, hop table check, CRC convention when one passes
            ok = [d for d in a1 if d.get("crc_ok")]
            out["fields"]["afhds1"] = {
                "packets": len(a1), "crc_ok": len(ok), "tx_ids": sorted({d["tx_id"] for d in (ok or a1)}),
                "hop_checks": sorted({d.get("hop_check") for d in a1 if d.get("hop_check")}),
                "hop_channels": a1[0]["hop_channels"], "binding": any(d["packet_name"] == "bind" for d in a1),
                "convention": ok[0]["convention"] if ok else None}
        eg = [e.decoded for e in em if e.label == "elrs_gfsk" and e.decoded]
        if eg:
            # ExpressLRS LR1121 GFSK (research, PLAN G3): verified SYNC packets -> UID bytes, OTA version, whitening
            out["fields"]["elrs_gfsk"] = {
                "sync_packets": len(eg), "uid4": eg[0]["uid4"], "uid5": eg[0]["uid5"],
                "sync_word": eg[0]["sync_word"], "ota": eg[0]["ota"], "ota_version": eg[0]["ota_version"],
                "whitening": sorted({d["whitening"] for d in eg}), "binding": any(d["binding"] for d in eg)}
        rl = [e.decoded for e in em if e.label == "radiolink" and e.decoded]
        if rl:
            # RadioLink (research, PLAN G4): sub-protocol, TX IDs, hop-table check, telemetry
            out["fields"]["radiolink"] = {
                "packets": len(rl), "sub_protocols": sorted({d["sub_protocol"] for d in rl}),
                "tx_ids": sorted({d["tx_id"] for d in rl}),
                "packet_kinds": sorted({d["packet_name"] for d in rl}),
                "hop_checks": sorted({d["hop_check"] for d in rl if d.get("hop_check")}),
                "hop_channels": next((d["hop_channels"] for d in rl if d.get("hop_channels")), None),
                "telemetry_requested": any(d.get("telemetry_request") for d in rl)}
        cc = [e.decoded for e in em if e.label in ("hott", "hitec", "sfhss") and e.decoded]
        if cc:
            # HoTT / Hitec / S-FHSS (CC2500; research, PLAN G12): TX IDs, packet kinds, hop table / index, model hint
            out["fields"]["cc2500_rc"] = {
                "packets": len(cc), "protocols": sorted({d["protocol"] for d in cc}),
                "tx_ids": sorted({d["tx_id"] for d in cc}),
                "rx_ids": sorted({d["rx_id"] for d in cc if d.get("rx_id")}),
                "packet_kinds": sorted({d["packet_name"] for d in cc}),
                "hop_table": next((d["hop_table"] for d in cc if d.get("hop_table") is not None), None),
                "hop_checks": sorted({d["hop_check"] for d in cc if d.get("hop_check")}),
                "model_hints": sorted({d["model_hint"] for d in cc if d.get("model_hint")})}
        hb = [e.decoded for e in em if e.label == "hubsan" and e.decoded]
        if hb:
            # Hubsan (research, PLAN G5): session ID codes, packet kinds, model from constant bytes
            out["fields"]["hubsan"] = {
                "packets": len(hb), "id_codes": sorted({d["id_code"] for d in hb}),
                "packet_kinds": sorted({d["packet_name"] for d in hb}),
                "channel_checks": sorted({d["channel_check"] for d in hb if d.get("channel_check")})}
        ty = [e.decoded for e in em if e.label == "toy_drone" and e.decoded]
        if ty:
            # toy-drone protocols on nRF24 / XN297 (research, PLAN G6): protocol, addresses, TX IDs
            out["fields"]["toy_drone"] = {
                "packets": len(ty), "protocols": sorted({d["protocol"] for d in ty}),
                "addresses": sorted({d["address"] for d in ty}), "tx_ids": sorted({d["tx_id"] for d in ty}),
                "checksum_ok": sum(bool(d.get("checksum_ok")) for d in ty)}
        xn = [e.decoded for e in em if e.label in ("bayang", "xn297_like") and e.decoded]
        if xn:
            # XN297 / Bayang (research, PLAN R8): addresses, Bayang TX IDs, packet kinds, CRC residuals
            ok = [d for d in xn if d.get("checksum_ok")]
            out["fields"]["xn297"] = {
                "packets": len(xn), "bayang_checksum_ok": len(ok), "addresses": sorted({d["address"] for d in xn}),
                "tx_ids": sorted({d["tx_id"] for d in ok if d.get("tx_id")}),
                "packet_kinds": sorted({d["packet_name"] for d in ok}),
                "crc_residuals": sorted({d["crc_residual"] for d in ok})}
        zb = [e.decoded for e in em if e.label.startswith("ieee802154") and e.decoded]
        if zb:
            # IEEE 802.15.4 (research, PLAN T4): CRC-verified frames, channel, PANs, addresses, vendors
            ok = [d for d in zb if d.get("crc_ok")]
            out["fields"]["ieee802154"] = {
                "frames": len(zb), "crc_ok": len(ok), "channels": sorted({d["channel"] for d in zb}),
                "frame_types": sorted({str(d["frame_type"]) for d in ok}),
                "pan_ids": sorted({d[k] for d in ok for k in ("dst_pan", "src_pan") if d.get(k)}),
                "addresses": sorted({d[k] for d in ok for k in ("dst_addr", "src_addr") if d.get(k)}),
                "vendors": sorted({d["vendor"] for d in ok if d.get("vendor")})}
        fl = [e.decoded for e in em if e.label == "flarm_like" and e.decoded]
        if fl:
            # FLARM (research, PLAN T5): clear header only; the payload is obfuscated and not decoded
            out["fields"]["flarm"] = {"beacons": len(fl),
                                      "radio_ids": sorted({d["radio_id"] for d in fl if d.get("radio_id")}),
                                      "address_types": sorted({d["address_type"] for d in fl if d.get("address_type")})}
        ml = [e for e in em if e.label.startswith("mlrs")]
        if ml:
            # mLRS (research, PLAN R3 / T6): air modes; CRC-verified frames give the per-binding FrameSyncWord
            ok = [e.decoded for e in ml if (e.decoded or {}).get("crc_ok")]
            out["fields"]["mlrs"] = {"packets": len(ml), "crc_ok": len(ok),
                                     "frame_sync_words": sorted({d["frame_sync_word"] for d in ok}),
                                     "modes": sorted({m for e in ml for m in e.features.get("mlrs_modes") or []})}
        frames = [e.features["frame"] for e in em if e.features.get("frame")]
        if frames:
            # blind FSK frame structure (research): sync-word candidates seen (a repeated one identifies the link),
            # known sync words named as "-like" candidates. Counted on the first 24 bits: the preamble run can absorb
            # the sync word's first bits, so the 32-bit window may end in payload bits
            from collections import Counter
            cnt = Counter(f["sync32_hex"][:6] for f in frames)
            known = sorted({f["known_sync"]["name"] for f in frames if f.get("known_sync")})
            out["fields"]["fsk_frame"] = {
                "packets": len(frames), "sync24_candidates": [[k, v] for k, v in cnt.most_common(4)],
                "preamble_bits": sorted({f["preamble_bits"] for f in frames}), "known_sync": known}
        if "model" not in out["fields"]:              # model level (research, PLAN M5)
            from .models import model
            mb = model(em)
            if mb:
                out["fields"]["model"] = mb
        coh = [e.decoded for e in em if (e.decoded or {}).get("crc_ok") and (e.decoded or {}).get("demod") == "coherent"]
        if coh:
            # low-SNR coherent GFSK path (research, PLAN L1): how many adverts, how many CRC Chase corrections, and
            # the measured modulation index (a transmitter property: BLE allows 0.45-0.55)
            import numpy as np
            out["fields"]["ble_coherent"] = {"adverts": len(coh), "bit_corrections": sum(d["bit_corrections"] for d in coh),
                                             "mod_index": round(float(np.median([d["mod_index"] for d in coh])), 2)}
        return out

    def modulation(self, em):
        m = em[0].features.get("modulation")
        return {"lora": "LoRa CSS"}.get(m, m or "GFSK")

    def frequency_label_check(self, signals):
        """A CRC-verified BLE advert passes only when de-whitened with its true channel
        index, so (labelled frequency - true advertising-channel frequency) is the
        capture's frequency-label error."""
        from .label_check import check_frequency_label
        return check_frequency_label([e for s in signals if s.features.get("decoder") == self.name
                                      for e in s.features.get("_emissions", [])])


PLUGINS = [RcLinkDecoder()]
