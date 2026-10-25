"""IEEE 802.11 MAC frame dissector (802.11-2020, 802.11ax/be element IDs, WPS 2.0, Wi-Fi P2P).

    parse(mpdu_without_fcs) -> dict     flat fields + nested element dicts + a one-line "summary"
    vendor(mac) -> str | None           IEEE registry name (MA-L/MA-M/MA-S); group / randomised flagged

Everything that a single receiver can read without keys is reported: frame control, address roles,
sequence, QoS; beacon/probe fixed fields (TSF -> AP uptime, interval, capabilities); every common element
(SSID, rates, DS channel, TIM, country, ERP, BSS load, RSN/WPA security, HT/VHT/HE/EHT capabilities and
operation, 802.11k/r/v/u hints, WMM, WPS device identity, Wi-Fi Direct device name); auth/assoc/deauth
codes; action categories; Block Ack; and for data frames the CCMP/TKIP/WEP header (PN, key ID), or when
unencrypted LLC/SNAP -> EAPOL handshake message / ARP / IPv4 / IPv6 headers.
"""
from __future__ import annotations

import gzip
from pathlib import Path

from spectrum_intel.core import mavlink

from . import droneid

FTYPE = {0: "mgmt", 1: "ctrl", 2: "data", 3: "ext"}
SUBTYPE = {
    (0, 0): "assoc-req", (0, 1): "assoc-resp", (0, 2): "reassoc-req", (0, 3): "reassoc-resp",
    (0, 4): "probe-req", (0, 5): "probe-resp", (0, 6): "timing-adv", (0, 8): "beacon", (0, 9): "atim",
    (0, 10): "disassoc", (0, 11): "auth", (0, 12): "deauth", (0, 13): "action", (0, 14): "action-no-ack",
    (1, 2): "trigger", (1, 4): "bf-report-poll", (1, 5): "ndp-announcement", (1, 6): "ctrl-frame-ext",
    (1, 7): "ctrl-wrapper", (1, 8): "block-ack-req", (1, 9): "block-ack", (1, 10): "ps-poll", (1, 11): "rts",
    (1, 12): "cts", (1, 13): "ack", (1, 14): "cf-end", (1, 15): "cf-end+cf-ack",
    (2, 0): "data", (2, 1): "data+cf-ack", (2, 2): "data+cf-poll", (2, 3): "data+cf-ack+cf-poll",
    (2, 4): "null", (2, 5): "cf-ack", (2, 6): "cf-poll", (2, 7): "cf-ack+cf-poll", (2, 8): "qos-data",
    (2, 9): "qos-data+cf-ack", (2, 10): "qos-data+cf-poll", (2, 11): "qos-data+cf-ack+cf-poll",
    (2, 12): "qos-null", (2, 14): "qos-cf-poll", (2, 15): "qos-cf-ack+cf-poll",
    (3, 0): "dmg-beacon", (3, 1): "s1g-beacon",
}
CAPS = [(0, "ESS"), (1, "IBSS"), (4, "privacy"), (5, "short-preamble"), (8, "spectrum-mgmt"), (9, "QoS"),
        (10, "short-slot"), (11, "APSD"), (12, "radio-measurement"), (13, "EPD"), (14, "delayed-BA"),
        (15, "immediate-BA")]
CIPHER = {0: "use-group", 1: "WEP-40", 2: "TKIP", 4: "CCMP-128", 5: "WEP-104", 6: "BIP-CMAC-128",
          7: "group-traffic-not-allowed", 8: "GCMP-128", 9: "GCMP-256", 10: "CCMP-256", 11: "BIP-GMAC-128",
          12: "BIP-GMAC-256", 13: "BIP-CMAC-256"}
AKM = {1: "802.1X", 2: "PSK", 3: "FT-802.1X", 4: "FT-PSK", 5: "802.1X-SHA256", 6: "PSK-SHA256", 7: "TDLS",
       8: "SAE", 9: "FT-SAE", 10: "AP-PeerKey", 11: "802.1X-SuiteB", 12: "802.1X-SuiteB-192",
       13: "FT-802.1X-SHA384", 14: "FILS-SHA256", 15: "FILS-SHA384", 16: "FT-FILS-SHA256",
       17: "FT-FILS-SHA384", 18: "OWE", 19: "FT-PSK-SHA384", 20: "PSK-SHA384", 22: "802.1X-SHA384",
       24: "SAE-EXT-KEY", 25: "FT-SAE-EXT-KEY"}
AUTH_ALG = {0: "open", 1: "shared-key", 2: "FT", 3: "SAE", 4: "FILS-SK", 5: "FILS-SK-PFS", 6: "FILS-PK",
            7: "PASN", 0xFFFF: "vendor"}
REASON = {1: "unspecified", 2: "previous authentication no longer valid", 3: "station leaving (deauthenticated)",
          4: "inactivity", 5: "AP cannot handle all associated stations", 6: "class 2 frame from non-authenticated STA",
          7: "class 3 frame from non-associated STA", 8: "station leaving BSS (disassociated)",
          9: "not authenticated", 10: "power capability unacceptable", 11: "supported channels unacceptable",
          12: "BSS transition management", 13: "invalid element", 14: "MIC failure", 15: "4-way handshake timeout",
          16: "group key handshake timeout", 17: "element in 4-way handshake differs", 18: "invalid group cipher",
          19: "invalid pairwise cipher", 20: "invalid AKMP", 21: "unsupported RSNE version",
          22: "invalid RSNE capabilities", 23: "802.1X authentication failed", 24: "cipher suite rejected",
          25: "TDLS teardown unreachable", 26: "TDLS teardown unspecified", 34: "excessive frame loss (poor channel)",
          36: "requesting STA leaving BSS", 39: "peer timeout", 45: "peer key mismatch", 66: "invalid PMKID"}
STATUS = {0: "success", 1: "unspecified failure", 10: "capabilities mismatch", 11: "reassociation denied",
          12: "association denied", 13: "authentication algorithm not supported", 14: "auth sequence out of order",
          15: "challenge failure", 16: "authentication timeout", 17: "AP full", 18: "basic rates mismatch",
          30: "rejected temporarily (try later)", 31: "robust management frame policy violation",
          37: "request declined", 40: "invalid element", 41: "invalid group cipher", 42: "invalid pairwise cipher",
          43: "invalid AKMP", 53: "invalid PMKID", 72: "invalid SAE group", 76: "anti-clogging token required",
          77: "finite cyclic group not supported", 82: "rejected with suggested BSS transition",
          123: "SAE hash-to-element required", 126: "SAE public key"}
ACTION_CAT = {0: "spectrum-mgmt", 1: "QoS", 2: "DLS", 3: "block-ack", 4: "public", 5: "radio-measurement",
              6: "fast-BSS-transition", 7: "HT", 8: "SA-query", 9: "protected-dual-of-public", 10: "WNM",
              11: "unprotected-WNM", 12: "TDLS", 13: "mesh", 14: "multihop", 15: "self-protected", 16: "DMG",
              18: "fast-session-transfer", 19: "robust-AV-streaming", 20: "unprotected-DMG", 21: "VHT",
              30: "HE", 31: "protected-HE", 33: "protected-FTM", 36: "EHT", 37: "protected-EHT",
              126: "vendor-protected", 127: "vendor"}
VENDOR_IE = {"00:50:f2": "Microsoft", "00:10:18": "Broadcom", "00:03:7f": "Qualcomm Atheros",
             "00:0c:43": "MediaTek (Ralink)", "00:0c:e7": "MediaTek", "00:e0:4c": "Realtek", "00:17:f2": "Apple",
             "50:6f:9a": "Wi-Fi Alliance", "00:90:4c": "Broadcom (Epigram)", "8c:fd:f0": "Qualcomm",
             "00:13:92": "Ruckus", "00:0b:86": "Aruba", "00:40:96": "Cisco", "00:a0:40": "Apple",
             "00:26:86": "Quantenna", "f4:f5:e8": "Google", "00:13:74": "Atheros", "00:15:6d": "Ubiquiti",
             "00:25:00": "Apple", "00:1d:0f": "TP-Link", "ac:85:3d": "Huawei", "00:e0:fc": "Huawei",
             "00:19:3b": "LG", "00:16:32": "Samsung", "00:00:f0": "Samsung"}
WFA_TYPES = {9: "Wi-Fi Direct (P2P)", 0x10: "Hotspot 2.0", 0x12: "OSEN", 0x16: "MBO/OCE",
             0x1C: "OWE transition", 0x1A: "Wi-Fi Easy Connect (DPP)", 0x1B: "Multi-AP (EasyMesh)"}
WPS_ATTR = {0x1011: "device_name", 0x1021: "manufacturer", 0x1023: "model_name", 0x1024: "model_number",
            0x1042: "serial_number"}
WPS_DEVICE_CAT = {1: "computer", 2: "input device", 3: "printer/scanner", 4: "camera", 5: "storage",
                  6: "network infrastructure (AP/router)", 7: "display", 8: "multimedia", 9: "gaming",
                  10: "telephone", 11: "audio", 12: "docking", 255: "other"}
ETHERTYPE = {0x0800: "IPv4", 0x0806: "ARP", 0x86DD: "IPv6", 0x888E: "EAPOL", 0x88C7: "802.11 pre-auth",
             0x8100: "VLAN", 0x88CC: "LLDP", 0x893A: "IEEE 1905 (EasyMesh)", 0x890D: "802.11 TDLS",
             0x88B4: "WAI (WAPI)"}
IP_PROTO = {1: "ICMP", 2: "IGMP", 6: "TCP", 17: "UDP", 58: "ICMPv6"}
PORTS = {53: "DNS", 67: "DHCP", 68: "DHCP", 123: "NTP", 137: "NetBIOS", 138: "NetBIOS", 1900: "SSDP",
         5353: "mDNS", 5355: "LLMNR", 443: "HTTPS", 80: "HTTP", 546: "DHCPv6", 547: "DHCPv6",
         14550: "MAVLink (GCS)", 14555: "MAVLink", 14540: "MAVLink (offboard)", 5760: "MAVLink (TCP)",
         8889: "Tello SDK command", 8890: "Tello SDK state", 11111: "Tello video"}
OFDM_RATES = {6.0, 9.0, 12.0, 18.0, 24.0, 36.0, 48.0, 54.0}

_OUI: dict | None = None


def _oui_db() -> dict:
    global _OUI
    if _OUI is None:
        _OUI = {}
        p = Path(__file__).with_name("data") / "oui.tsv.gz"
        try:
            with gzip.open(p, "rt", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    k, _, v = line.rstrip("\n").partition("\t")
                    _OUI[k] = v
        except OSError:
            pass
    return _OUI


def vendor(mac: str | None) -> str | None:
    """Registry vendor of a unicast, universally administered address; else what the address is."""
    if not mac:
        return None
    o1 = int(mac[:2], 16)
    if o1 & 0x01:                                                   # I/G bit: group address, not a device
        return "broadcast" if mac.lower() == "ff:ff:ff:ff:ff:ff" else "multicast"
    if o1 & 0x02:
        return "locally-administered (randomized/private)"
    h = mac.replace(":", "").upper()
    db = _oui_db()
    for n in (9, 7, 6):
        if h[:n] in db:
            return db[h[:n]]
    return "unknown OUI"


def oui_name(oui: str) -> str:
    """Name of a 3-byte organisation identifier (vendor-specific elements)."""
    o = oui.lower()
    if o in VENDOR_IE:
        return VENDOR_IE[o]
    return _oui_db().get(o.replace(":", "").upper(), o)


def _mac(b):
    return ":".join(f"{x:02x}" for x in b)


def _u16(b, i=0):
    return b[i] | (b[i + 1] << 8)


def _uptime(us):
    s = us // 1_000_000
    d, s = divmod(s, 86400)
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    return f"{d}d {h:02d}:{m:02d}:{s:02d}"


# ------------------------------------------------------------------------------------------- elements
def _suite(b, i):
    return _mac(b[i:i + 3]), b[i + 3]


def _rsn(v, oui_expect="00:0f:ac", wpa1=False):
    """RSN (or Microsoft WPA1) element body -> dict."""
    out = {"version": _u16(v, 0)}
    i = 2
    try:
        if i + 4 <= len(v):
            o, t = _suite(v, i)
            out["group_cipher"] = CIPHER.get(t, f"{o}:{t}") if o == oui_expect else f"{o}:{t}"
            i += 4
        if i + 2 <= len(v):
            n = _u16(v, i)
            i += 2
            out["pairwise_ciphers"] = [CIPHER.get(_suite(v, i + 4 * k)[1], "?") for k in range(n)]
            i += 4 * n
        if i + 2 <= len(v):
            n = _u16(v, i)
            i += 2
            akm = []
            for k in range(n):
                o, t = _suite(v, i + 4 * k)
                akm.append(AKM.get(t, f"{t}") if o == oui_expect else f"{oui_name(o)}:{t}")
            out["akm"] = akm
            i += 4 * n
        if not wpa1 and i + 2 <= len(v):
            cap = _u16(v, i)
            out["mfp_required"] = bool(cap & 0x40)
            out["mfp_capable"] = bool(cap & 0x80)
            i += 2
            if i + 2 <= len(v):
                n = _u16(v, i)
                out["n_pmkid"] = n
                i += 2 + 16 * n
            if i + 4 <= len(v):
                out["group_mgmt_cipher"] = CIPHER.get(v[i + 3], "?")
    except IndexError:
        out["truncated"] = True
    return out


def _wps(v):
    out = {}
    i = 0
    while i + 4 <= len(v):
        t, n = (v[i] << 8) | v[i + 1], (v[i + 2] << 8) | v[i + 3]
        val = v[i + 4:i + 4 + n]
        if t in WPS_ATTR:
            out[WPS_ATTR[t]] = val.decode("utf-8", "replace").strip("\x00 ")
        elif t == 0x1044 and n >= 1:
            out["state"] = {1: "not configured", 2: "configured"}.get(val[0], val[0])
        elif t == 0x1047 and n == 16:
            out["uuid"] = val.hex()
        elif t == 0x1054 and n >= 8:
            out["device_category"] = WPS_DEVICE_CAT.get((val[0] << 8) | val[1], (val[0] << 8) | val[1])
        elif t == 0x103B and n >= 1:
            out["response_type"] = {0: "enrollee (open 802.1X)", 1: "enrollee", 2: "registrar", 3: "AP"}.get(val[0])
        elif t == 0x1041 and n >= 1:
            out["selected_registrar"] = bool(val[0])
        elif t == 0x1012 and n >= 2:
            out["device_password_id"] = {0: "PIN", 4: "push-button"}.get((val[0] << 8) | val[1], (val[0] << 8) | val[1])
        i += 4 + n
    return out


def _p2p(v):
    """Wi-Fi Direct attributes: device address and name (attribute 13, Device Info)."""
    out = {}
    i = 0
    while i + 3 <= len(v):
        a, n = v[i], _u16(v, i + 1)
        val = v[i + 3:i + 3 + n]
        if a == 13 and n >= 17:
            out["device_address"] = _mac(val[0:6])
            nsec = val[16]
            j = 17 + 8 * nsec
            if j + 4 <= len(val) and ((val[j] << 8) | val[j + 1]) == 0x1011:
                ln = (val[j + 2] << 8) | val[j + 3]
                out["device_name"] = val[j + 4:j + 4 + ln].decode("utf-8", "replace")
            cat = (val[8] << 8) | val[9]
            out["device_category"] = WPS_DEVICE_CAT.get(cat, cat)
        elif a == 3 and n >= 6:
            out["device_id"] = _mac(val[0:6])
        i += 3 + n
    return out


def _rates(v):
    return [((b & 0x7F) / 2, bool(b & 0x80)) for b in v if (b & 0x7F) not in (0x7F, 0x7E)]


def parse_elements(tlv: bytes) -> dict:
    """Tagged elements -> dict (unknown IDs are listed in element_ids)."""
    out: dict = {"element_ids": []}
    rates, basic, vend, ext_ids = [], [], [], []
    i = 0
    while i + 2 <= len(tlv):
        tag, ln = tlv[i], tlv[i + 1]
        v = tlv[i + 2:i + 2 + ln]
        if len(v) < ln:
            out["elements_truncated"] = True
            break
        out["element_ids"].append(tag)
        try:
            if tag == 0:
                if ln == 0 or all(b == 0 for b in v):
                    out["ssid"] = "<broadcast>" if ln == 0 else "<hidden>"
                    out["ssid_hidden"] = True
                else:
                    out["ssid"] = v.decode("utf-8", "replace")
            elif tag in (1, 50):
                for r, b in _rates(v):
                    rates.append(r)
                    if b:
                        basic.append(r)
            elif tag == 3 and ln >= 1:
                out["ds_channel"] = int(v[0])
            elif tag == 5 and ln >= 4:
                bc, pvb = v[2], v[3:]
                offset = (bc >> 1) * 2
                aids = [8 * (offset + k) + j for k, b in enumerate(pvb) for j in range(8) if b >> j & 1]
                aids = [a for a in aids if a]
                out["tim"] = {"dtim_count": v[0], "dtim_period": v[1], "multicast_buffered": bool(bc & 1),
                              "aids_with_traffic": aids[:32]}
            elif tag == 7 and ln >= 3:
                trip = [(v[k], v[k + 1], v[k + 2] if v[k + 2] < 128 else v[k + 2] - 256)
                        for k in range(3, ln - 2, 3) if v[k] < 201]
                out["country"] = {"code": v[:2].decode("ascii", "replace"),
                                  "environment": {0x20: "any", 0x4F: "outdoor", 0x49: "indoor",
                                                  0x58: "non-country"}.get(v[2], hex(v[2])),
                                  "channels": [{"first": a, "n": b, "max_dbm": c} for a, b, c in trip]}
            elif tag == 11 and ln >= 5:
                out["bss_load"] = {"stations": _u16(v, 0), "channel_utilisation_pct": round(v[2] / 255 * 100, 1),
                                   "available_admission_capacity": _u16(v, 3)}
            elif tag == 32 and ln >= 1:
                out["power_constraint_db"] = v[0]
            elif tag == 42 and ln >= 1:
                out["erp"] = {"non_erp_present": bool(v[0] & 1), "use_protection": bool(v[0] & 2),
                              "barker_preamble_mode": bool(v[0] & 4)}
            elif tag == 45 and ln >= 26:
                cap = _u16(v, 0)
                mcs = v[3:7]
                out["ht"] = {"40mhz": bool(cap & 0x02), "sgi20": bool(cap & 0x20), "sgi40": bool(cap & 0x40),
                             "greenfield": bool(cap & 0x10), "ldpc": bool(cap & 0x01),
                             "spatial_streams": sum(1 for b in mcs if b)}
            elif tag == 61 and ln >= 2:
                sec = {0: "none", 1: "above", 3: "below"}.get(v[1] & 3, "?")
                out["ht_operation"] = {"primary_channel": v[0], "secondary_channel": sec,
                                       "width_mhz": 40 if (v[1] & 4 and sec != "none") else 20}
            elif tag == 48:
                out["rsn"] = _rsn(v)
            elif tag == 54 and ln >= 2:
                out["mobility_domain"] = f"{_u16(v, 0):04x}"
            elif tag == 70:
                out["radio_measurement_802_11k"] = True
            elif tag == 107:
                out["interworking_802_11u"] = True
            elif tag == 127:
                ext = int.from_bytes(v, "little")
                out["ext_capabilities"] = [n for b, n in ((19, "BSS-transition (802.11v)"), (31, "interworking"),
                                                          (62, "opmode-notification"), (70, "FTM-responder"),
                                                          (71, "FTM-initiator"), (77, "TWT-requester"),
                                                          (78, "TWT-responder")) if ext >> b & 1]
            elif tag == 191 and ln >= 12:
                cap = int.from_bytes(v[0:4], "little")
                rx_map = _u16(v, 4)
                ss = sum(1 for k in range(8) if (rx_map >> (2 * k)) & 3 != 3)
                out["vht"] = {"width_set": {0: "80", 1: "160", 2: "160/80+80"}.get((cap >> 2) & 3, "?"),
                              "sgi80": bool(cap & 0x20), "spatial_streams": ss,
                              "su_beamformer": bool(cap & (1 << 11)), "mu_beamformer": bool(cap & (1 << 19))}
            elif tag == 192 and ln >= 3:
                out["vht_operation"] = {"width": {0: "20/40", 1: "80", 2: "160", 3: "80+80"}.get(v[0], v[0]),
                                        "center_seg0": v[1], "center_seg1": v[2]}
            elif tag == 221 and ln >= 4:
                oui, typ = _mac(v[0:3]), v[3]
                name = oui_name(oui)
                if oui == "00:50:f2" and typ == 1:
                    out["wpa"] = _rsn(v[4:], "00:50:f2", wpa1=True)
                    name = "Microsoft WPA"
                elif oui == "00:50:f2" and typ == 2:
                    out["wmm"] = True
                    name = "Microsoft WMM/WME"
                elif oui == "00:50:f2" and typ == 4:
                    out["wps"] = _wps(v[4:])
                    name = "Microsoft WPS"
                elif oui == "50:6f:9a":
                    name = f"Wi-Fi Alliance {WFA_TYPES.get(typ, typ)}"
                    if typ == 9:
                        out["p2p"] = _p2p(v[4:])
                    elif typ == 0x1C and ln >= 11:
                        out["owe_transition_bssid"] = _mac(v[4:10])
                    elif typ == 0x10:
                        out["hotspot2"] = True
                else:
                    rid = droneid.from_vendor_ie(oui, typ, v[4:])     # Remote ID (F3411 / DJI)
                    if rid:
                        out["remote_id"] = rid
                        name = rid["standard"]
                vend.append(name)
            elif tag == 255 and ln >= 1:
                ext_ids.append(v[0])
                if v[0] == 35:
                    out["he"] = True
                elif v[0] == 36 and ln >= 7:
                    out["he_operation"] = {"bss_color": v[4] & 0x3F}
                elif v[0] in (106, 108):
                    out["eht"] = True
        except (IndexError, ValueError):
            out.setdefault("malformed_elements", []).append(tag)
        i += 2 + ln
    if rates:
        out["rates_mbps"] = sorted(set(rates))
        out["basic_rates_mbps"] = sorted(set(basic))
    if vend:
        out["vendor_elements"] = vend
    if ext_ids:
        out["ext_element_ids"] = ext_ids
    return out


def security(info: dict) -> str | None:
    """Human summary from RSN/WPA elements and the privacy bit (beacons, probe responses, assoc frames)."""
    rsn, wpa = info.get("rsn"), info.get("wpa")
    priv = "privacy" in info.get("capabilities", [])
    if not rsn and not wpa:
        if "capabilities" not in info:
            return None
        return "WEP" if priv else "open"
    parts = []
    if rsn:
        akm = set(rsn.get("akm", []))
        ci = "/".join(rsn.get("pairwise_ciphers", [])) or rsn.get("group_cipher", "?")
        if "OWE" in akm:
            kind = "OWE (Enhanced Open)"
        elif akm & {"SAE", "FT-SAE", "SAE-EXT-KEY", "FT-SAE-EXT-KEY"}:
            kind = "WPA2/WPA3-Personal (transition)" if akm & {"PSK", "FT-PSK", "PSK-SHA256"} else "WPA3-Personal (SAE)"
        elif akm & {"802.1X-SuiteB-192", "802.1X-SHA384"}:
            kind = "WPA3-Enterprise 192-bit"
        elif any("802.1X" in a or "FILS" in a for a in akm):
            kind = "WPA2-Enterprise (802.1X)"
        elif akm & {"PSK", "FT-PSK", "PSK-SHA256"}:
            kind = "WPA2-Personal (PSK)"
        else:
            kind = "RSN " + ",".join(sorted(akm))
        parts.append(f"{kind} {ci}" + (" +FT" if any(a.startswith("FT") for a in akm) else "")
                     + (" MFP-required" if rsn.get("mfp_required") else ""))
    if wpa:
        parts.append("WPA1 " + "/".join(wpa.get("pairwise_ciphers", [])))
    return "; ".join(parts)


def phy_generation(info: dict) -> str | None:
    rates = set(info.get("rates_mbps", []))
    if not rates and not info.get("ht"):
        return None
    std = ["b"] if rates & {1.0, 2.0, 5.5, 11.0} else []
    if rates & OFDM_RATES:
        std.append("g")
    for k, n in (("ht", "n"), ("vht", "ac"), ("he", "ax"), ("eht", "be")):
        if info.get(k):
            std.append(n)
    return "802.11" + "/".join(std) if std else None


# ------------------------------------------------------------------------------------------- data bodies
def _llc(body: bytes) -> dict:
    out = {}
    if len(body) < 8 or body[0:3] != b"\xaa\xaa\x03":
        return out
    et = (body[6] << 8) | body[7]
    out["ethertype"] = f"0x{et:04x}"
    out["protocol"] = ETHERTYPE.get(et, out["ethertype"])
    p = body[8:]
    try:
        if et == 0x888E and len(p) >= 4:
            ev = {"eapol_version": p[0], "eapol_type": {0: "EAP", 1: "start", 2: "logoff", 3: "key"}.get(p[1], p[1])}
            if p[1] == 3 and len(p) >= 99:
                ki = (p[5] << 8) | p[6]
                pair, inst, ack, mic, sec = ki & 8, ki & 0x40, ki & 0x80, ki & 0x100, ki & 0x200
                if pair:
                    msg = 1 if ack and not mic else 3 if ack and mic and inst else \
                        4 if mic and sec and not ack else 2 if mic and not ack else None
                    ev["handshake"] = f"4-way handshake message {msg}" if msg else "pairwise key"
                else:
                    ev["handshake"] = "group key handshake"
                ev["key_descriptor_version"] = ki & 7
                ev["replay_counter"] = int.from_bytes(p[9:17], "big")
            elif p[1] == 0 and len(p) >= 9:
                code = {1: "request", 2: "response", 3: "success", 4: "failure"}.get(p[4], p[4])
                ev["eap"] = code
                if p[4] in (1, 2) and len(p) >= 9:
                    ev["eap_method"] = {1: "identity", 13: "EAP-TLS", 21: "EAP-TTLS", 25: "PEAP",
                                        43: "EAP-FAST", 52: "EAP-PWD"}.get(p[8], p[8])
                    if p[8] == 1 and p[4] == 2:
                        ev["identity"] = p[9:4 + ((p[6] << 8) | p[7])].decode("utf-8", "replace")
            out["eapol"] = ev
        elif et == 0x0806 and len(p) >= 28:
            out["arp"] = {"op": {1: "request", 2: "reply"}.get(_be16(p, 6), _be16(p, 6)),
                          "sender_mac": _mac(p[8:14]), "sender_ip": _ip4(p[14:18]),
                          "target_mac": _mac(p[18:24]), "target_ip": _ip4(p[24:28])}
        elif et == 0x0800 and len(p) >= 20:
            ihl = (p[0] & 15) * 4
            proto = p[9]
            ip = {"src": _ip4(p[12:16]), "dst": _ip4(p[16:20]), "proto": IP_PROTO.get(proto, proto), "ttl": p[8]}
            if proto in (6, 17) and len(p) >= ihl + 4:
                sp, dp = _be16(p, ihl), _be16(p, ihl + 2)
                ip.update(sport=sp, dport=dp)
                svc = PORTS.get(dp) or PORTS.get(sp)
                if svc:
                    ip["service"] = svc
                # MAVLink telemetry over an open Wi-Fi link (UDP 14550 bridges, ESP8266 / ESP32 DroneBridge in open
                # mode; research, PLAN G2): any CRC-valid frame in the UDP / TCP payload (CRC_EXTRA verifies it)
                off = ihl + 8 if proto == 17 else ihl + 4 * (p[ihl + 12] >> 4) if len(p) > ihl + 12 else len(p)
                mv = _mavlink(p[off:])
                if mv:
                    ip["mavlink"] = mv
                if proto == 17:                              # PLAN W4: Tello SDK, video streams (open links)
                    ip.update(_udp_app(p[off:], sp, dp))
            out["ipv4"] = ip
        elif et == 0x86DD and len(p) >= 40:
            nh = p[6]
            ip = {"src": _ip6(p[8:24]), "dst": _ip6(p[24:40]), "next_header": IP_PROTO.get(nh, nh)}
            if nh in (6, 17) and len(p) >= 44:
                ip.update(sport=_be16(p, 40), dport=_be16(p, 42))
            out["ipv6"] = ip
    except IndexError:
        out["truncated"] = True
    return out


def _be16(b, i):
    return (b[i] << 8) | b[i + 1]


def _ip4(b):
    return ".".join(str(x) for x in b)


def _ip6(b):
    return ":".join(f"{_be16(b, i):x}" for i in range(0, 16, 2))


def _crypto_header(body: bytes) -> dict:
    if len(body) < 4:
        return {}
    ext_iv = bool(body[3] & 0x20)
    key_id = body[3] >> 6
    if not ext_iv:
        return {"cipher": "WEP", "key_id": key_id, "iv": body[0:3].hex()}
    if len(body) < 8:
        return {"key_id": key_id}
    if body[1] == ((body[0] | 0x20) & 0x7F):                           # TKIP WEPSeed[1]
        pn = body[2] | (body[0] << 8) | (int.from_bytes(body[4:8], "little") << 16)
        return {"cipher": "TKIP", "key_id": key_id, "pn": pn}
    pn = body[0] | (body[1] << 8) | (int.from_bytes(body[4:8], "little") << 16)
    return {"cipher": "CCMP/GCMP", "key_id": key_id, "pn": pn}


# ------------------------------------------------------------------------------------------- frame
# wfb-ng / OpenIPC / OpenHD raw-802.11 links (research, PLAN G1): broadcast data frames, TA = BSSID = 57:42 + the
# 32-bit channel_id (link_id << 8 | radio_port, big-endian); payload type 0x01 data (8-byte nonce = block << 8 | frag)
# or 0x02 session (24-byte nonce + sealed session key, about once per second). link_id 0x7505D6 = the default
# link_domain "default" (wfb-ng and OpenIPC images). Content is ChaCha20-Poly1305 encrypted (not decrypted here).
WFB_DEFAULT_LINK_ID = 0x7505D6
WFB_TYPES = {1: "data", 2: "session"}


def _wfb_stream(port: int) -> str:
    base = port & 0x7F
    s = "video / data" if base < 16 else "MAVLink" if base < 32 else "tunnel" if base < 48 else \
        "bind" if base == 0x7F else "other"
    return s + (" (reverse: ground -> air)" if port & 0x80 else "")


def _wfb(a2: bytes, body: bytes) -> dict:
    link, port = int.from_bytes(a2[2:5], "big"), a2[5]
    out = {"link_id": link, "link_id_hex": f"0x{link:06X}", "default_link_domain": link == WFB_DEFAULT_LINK_ID,
           "radio_port": port, "stream": _wfb_stream(port)}
    if body:
        out["packet_type"] = WFB_TYPES.get(body[0], f"0x{body[0]:02x}")
        if body[0] == 1 and len(body) >= 9:
            n = int.from_bytes(body[1:9], "big")
            out.update(block_idx=n >> 8, fragment_idx=n & 0xFF)
    return out


# ESP-NOW (Espressif vendor action frames) and the ExpressLRS Backpack inside them (research, PLAN G7): category 127,
# OUI 18:fe:34, 4 random bytes, element 221 (OUI 18:fe:34, type 4, version), body = MSP v2 for the Backpack
# ('$' 'X' '<' / '>' flags function(LE16) size(LE16) payload crc8-DVB-S2). The Backpack's MAC (sender and peer) is the
# ELRS bind UID with byte 0's multicast bit cleared.
MSP_NAMES = {0x0009: "ELRS bind", 0x000B: "ELRS request VTX packet", 0x0010: "backpack version",
             0x0011: "backpack CRSF telemetry", 0x0058: "VTX config (get)", 0x0059: "VTX config (set)",
             0x00B5: "ELRS set send UID", 0x00B6: "ELRS set OSD", 0x0300: "get channel index",
             0x0301: "set channel index", 0x0302: "get frequency", 0x0303: "set frequency",
             0x0304: "get recording state", 0x0305: "set recording state", 0x0306: "get VRX mode",
             0x0307: "set VRX mode", 0x0308: "get RSSI", 0x0309: "get battery voltage", 0x030A: "get firmware",
             0x030D: "set head tracking", 0x030E: "set RTC", 0x0380: "set mode", 0x0381: "get version",
             0x0382: "get status", 0x0383: "set pan / tilt / roll"}
VTX_BANDS = "ABEFRL"


def crc8_dvb_s2(data: bytes) -> int:
    c = 0
    for b in data:
        c ^= b
        for _ in range(8):
            c = ((c << 1) ^ 0xD5) & 0xFF if c & 0x80 else (c << 1) & 0xFF
    return c


def _msp(b: bytes):
    if len(b) < 9 or b[0:2] != b"$X" or b[2] not in b"<>":
        return None
    func, size = b[4] | b[5] << 8, b[6] | b[7] << 8
    if len(b) < 9 + size:
        return None
    pay = b[8:8 + size]
    out = {"direction": "command" if b[2] == 0x3C else "response", "function": f"0x{func:04X}",
           "function_name": MSP_NAMES.get(func, "?"), "payload_hex": pay.hex(),
           "crc_ok": crc8_dvb_s2(b[3:8 + size]) == b[8 + size]}
    idx = None
    if func == 0x0059 and len(pay) >= 2:
        v = pay[0] | pay[1] << 8
        if v < 48:
            idx = v
        else:
            out["vtx_frequency_mhz"] = v
    elif func == 0x0301 and len(pay) >= 1:
        idx = pay[0]
    if idx is not None and idx < 48:
        out.update(vtx_band=VTX_BANDS[idx // 8], vtx_channel=idx % 8 + 1)
    return out


def _espnow(b: bytes, info: dict):
    """b = the action frame body (category byte first)."""
    out = {"random": b[4:8].hex()}
    i = 8
    while i + 2 <= len(b):
        eid, ln = b[i], b[i + 1]
        v = b[i + 2:i + 2 + ln]
        if eid == 221 and len(v) >= 5 and _mac(v[0:3]) == "18:fe:34" and v[3] == 4:
            out.update(version=v[4], body_len=len(v) - 5)
            m = _msp(v[5:])
            if m:
                out["msp"] = m
            break
        i += 2 + ln
    m = out.get("msp")
    a1, a2 = info.get("addr1"), info.get("addr2")
    if m and m["crc_ok"] and a2 and (a1 == a2 or (a1 == "ff:ff:ff:ff:ff:ff" and m["function"] == "0x0009")):
        info["elrs_backpack"] = {"bind_uid": a2, "msp": m["function_name"],
                                 **{k: m[k] for k in ("vtx_band", "vtx_channel", "vtx_frequency_mhz") if k in m}}
    return out


# Tello SDK (Ryze / DJI Tello SDK 2.0 / 3.0 user guide): commands to UDP 8889 on the drone (192.168.10.1), its
# state string to UDP 8890, raw H.264 video to UDP 11111. Tello's own network is open by default.
TELLO_CMDS = {"command", "takeoff", "land", "streamon", "streamoff", "emergency", "up", "down", "left", "right",
              "forward", "back", "cw", "ccw", "flip", "go", "stop", "curve", "jump", "speed", "rc", "wifi", "mon",
              "moff", "mdirection", "ap", "motoron", "motoroff", "throwfly", "reboot", "setfps", "setbitrate",
              "setresolution", "port", "multiwifi", "downvision", "battery?", "speed?", "time?", "wifi?", "sdk?",
              "sn?", "hardware?", "active?", "ap?", "ssid?"}
TELLO_STATE = {"pitch": "pitch_deg", "roll": "roll_deg", "yaw": "yaw_deg", "vgx": "vgx_dm_s", "vgy": "vgy_dm_s",
               "vgz": "vgz_dm_s", "templ": "temp_low_c", "temph": "temp_high_c", "tof": "tof_cm", "h": "height_cm",
               "bat": "battery_pct", "baro": "baro_m", "time": "motor_time_s", "agx": "agx", "agy": "agy",
               "agz": "agz", "mid": "mission_pad"}
H264_NAL = {1: "slice", 5: "IDR slice", 6: "SEI", 7: "SPS", 8: "PPS", 9: "AUD"}
H265_NAL = {1: "slice", 19: "IDR slice", 20: "IDR slice", 32: "VPS", 33: "SPS", 34: "PPS", 35: "AUD", 39: "SEI"}


def _nal(b: bytes):
    """Annex-B start code at offset 0 + a valid H.264 / H.265 NAL header -> (codec, nal type name), else None."""
    sc = 4 if b[:4] == b"\0\0\0\1" else 3 if b[:3] == b"\0\0\1" else 0
    if not sc or len(b) < sc + 2 or b[sc] & 0x80:
        return None
    h = b[sc]
    if not h & 1 and (h >> 1) & 0x3F in H265_NAL and b[sc + 1] in (1, 2, 3):   # H.265: layer 0, TID 1..3
        return "H.265", H265_NAL[(h >> 1) & 0x3F]
    if h & 0x1F in H264_NAL and (h & 0x1F not in (7, 8, 5) or h & 0x60):  # SPS / PPS / IDR carry nal_ref_idc > 0
        return "H.264", H264_NAL[h & 0x1F]
    return None


def _udp_app(d: bytes, sp: int, dp: int) -> dict:
    """Application payloads of open drone links (research, PLAN W4): Tello SDK command / state, and video streams
    (raw H.264 / H.265 Annex B, or RTP with a dynamic payload type carrying H.264 NAL units)."""
    out = {}
    if 8889 in (sp, dp) and 0 < len(d) <= 128:
        try:
            t = d.decode("ascii").strip()
        except UnicodeDecodeError:
            t = ""
        w = t.split(" ")[0].lower() if t and t.isprintable() else ""
        if w in TELLO_CMDS or (sp == 8889 and (w in ("ok", "error") or w.replace(".", "", 1).isdigit())):
            if w in ("wifi", "ap", "multiwifi"):                         # carries an SSID and a password
                t = w + " <credentials not reported>"
            out["tello"] = {"port": "command", "direction": "to drone" if dp == 8889 else "from drone", "text": t}
    if 8890 in (sp, dp) and len(d) <= 300:
        kv = {}
        for item in d.decode("ascii", "replace").strip().split(";"):
            k, _, v = item.partition(":")
            if k in TELLO_STATE and v:
                try:
                    kv[TELLO_STATE[k]] = float(v) if "." in v else int(v)
                except ValueError:
                    pass
        if len(kv) >= 8 and {"battery_pct", "height_cm", "pitch_deg"} <= set(kv):
            out["tello"] = {"port": "state", "direction": "from drone", **kv}
    nal = _nal(d)
    if nal:
        out["video"] = {"transport": "raw Annex B", "codec": nal[0], "nal": nal[1]}
    elif len(d) >= 14 and d[0] >> 6 == 2 and 96 <= (d[1] & 0x7F) <= 127:
        cc = d[0] & 0x0F
        off = 12 + 4 * cc
        if len(d) > off + 1 and not d[off] & 0x80:
            t = d[off] & 0x1F
            if 1 <= t <= 23 or t in (24, 28):
                out["video"] = {"transport": "RTP", "codec": "H.264-like", "payload_type": d[1] & 0x7F,
                                "nal": H264_NAL.get(t, {24: "STAP-A", 28: "FU-A"}.get(t, f"type {t}"))}
    return out


def _mavlink(data: bytes):
    """CRC-valid MAVLink frames in a datagram -> {frames, summary} or None."""
    if len(data) < 8:
        return None
    fr = [f for f in mavlink.scan(data) if f.get("crc_ok")]
    return {"frames": len(fr), **mavlink.summary(fr)} if fr else None


def parse(mpdu) -> dict:
    """Dissect one MPDU (without FCS)."""
    mpdu = bytes(bytearray(int(v) & 0xFF for v in mpdu))
    if len(mpdu) < 2:
        return {"parse": "too short"}
    fc0, fc1 = mpdu[0], mpdu[1]
    ftype, sub = (fc0 >> 2) & 3, (fc0 >> 4) & 0xF
    to_ds, from_ds = fc1 & 1, (fc1 >> 1) & 1
    info = {"type": FTYPE.get(ftype, "?"), "subtype": SUBTYPE.get((ftype, sub), f"{ftype}/{sub}"),
            "protocol_version": fc0 & 3, "to_ds": to_ds, "from_ds": from_ds, "more_frag": (fc1 >> 2) & 1,
            "retry": (fc1 >> 3) & 1, "pwr_mgmt": (fc1 >> 4) & 1, "more_data": (fc1 >> 5) & 1,
            "protected": (fc1 >> 6) & 1, "order": (fc1 >> 7) & 1, "length_bytes": len(mpdu) + 4}
    if len(mpdu) >= 4:
        d = _u16(mpdu, 2)
        if info["subtype"] == "ps-poll":
            info["aid"] = d & 0x3FFF
        else:
            info["duration_us"] = d if d < 0x8000 else None
    st = info["subtype"]
    if ftype == 1:                                                  # control
        if len(mpdu) >= 10:
            info["addr1_ra"] = info["ra"] = _mac(mpdu[4:10])
        if st in ("rts", "block-ack-req", "block-ack", "ps-poll", "bf-report-poll", "ndp-announcement",
                  "cf-end", "trigger") and len(mpdu) >= 16:
            info["addr2_ta"] = info["ta"] = _mac(mpdu[10:16])
        if st == "ps-poll":
            info["bssid"] = info.get("ra")
        if st in ("block-ack-req", "block-ack") and len(mpdu) >= 20:
            ctl = _u16(mpdu, 16)
            info["block_ack"] = {"tid": ctl >> 12, "start_seq": _u16(mpdu, 18) >> 4,
                                 "bitmap": mpdu[20:28].hex() if st == "block-ack" else None}
        _vendors(info)
        info["summary"] = summary(info)
        return info
    if len(mpdu) >= 10:
        info["addr1"] = _mac(mpdu[4:10])
    if len(mpdu) >= 16:
        info["addr2"] = _mac(mpdu[10:16])
    if len(mpdu) >= 22:
        info["addr3"] = _mac(mpdu[16:22])
    if len(mpdu) >= 24:
        sc = _u16(mpdu, 22)
        info["seq"] = sc >> 4
        info["frag"] = sc & 0xF
    hdr = 24
    if ftype == 2 and to_ds and from_ds and len(mpdu) >= 30:
        info["addr4"] = _mac(mpdu[24:30])
        hdr = 30
    info["ra"], info["ta"] = info.get("addr1"), info.get("addr2")
    if ftype == 2:
        a1, a2, a3, a4 = (info.get(k) for k in ("addr1", "addr2", "addr3", "addr4"))
        if not to_ds and not from_ds:
            info.update(da=a1, sa=a2, bssid=a3, direction="IBSS/direct")
        elif not to_ds and from_ds:
            info.update(da=a1, bssid=a2, sa=a3, direction="AP -> STA")
        elif to_ds and not from_ds:
            info.update(bssid=a1, sa=a2, da=a3, direction="STA -> AP")
        else:
            info.update(da=a3, sa=a4, direction="WDS / mesh (4-address)")
        if sub & 8 and len(mpdu) >= hdr + 2:                         # QoS control
            q = _u16(mpdu, hdr)
            info["qos"] = {"tid": q & 15, "ack_policy": {0: "normal", 1: "no-ack", 2: "no-explicit", 3: "block-ack"}
                           [(q >> 5) & 3], "amsdu": bool(q & 0x80)}
            hdr += 2
        if info["order"] and len(mpdu) >= hdr + 4:
            hdr += 4                                                # HT control
        body = mpdu[hdr:]
        info["body_bytes"] = len(body)
        if (sub == 0 and info.get("addr1") == "ff:ff:ff:ff:ff:ff" and str(info.get("addr2", "")).startswith("57:42:")
                and info.get("addr3") == info.get("addr2")):
            info["wfb_ng"] = _wfb(mpdu[10:16], body)
        if info["protected"]:
            info["encryption"] = _crypto_header(body)
        elif body and not (sub & 4):                                 # not a null frame
            if info.get("qos", {}).get("amsdu"):
                info["amsdu"] = True
            else:
                info.update(_llc(body))
    elif ftype == 0:
        info["bssid"] = info.get("addr3")
        if info["protected"]:
            info["encryption"] = _crypto_header(mpdu[24:])
            _vendors(info)
            info["summary"] = summary(info)
            return info
        b = mpdu[24:]
        ie = None
        if st in ("beacon", "probe-resp") and len(b) >= 12:
            tsf = int.from_bytes(b[0:8], "little")
            info.update(tsf_us=tsf, ap_uptime=_uptime(tsf), beacon_interval_tu=_u16(b, 8),
                        beacon_interval_ms=round(_u16(b, 8) * 1.024, 1))
            info["capabilities"] = _caps(_u16(b, 10))
            ie = b[12:]
        elif st == "probe-req":
            ie = b
        elif st in ("assoc-req", "reassoc-req") and len(b) >= 4:
            info["capabilities"] = _caps(_u16(b, 0))
            info["listen_interval"] = _u16(b, 2)
            if st == "reassoc-req" and len(b) >= 10:
                info["current_ap"] = _mac(b[4:10])
                ie = b[10:]
            else:
                ie = b[4:]
        elif st in ("assoc-resp", "reassoc-resp") and len(b) >= 6:
            info["capabilities"] = _caps(_u16(b, 0))
            info["status"], info["status_text"] = _u16(b, 2), STATUS.get(_u16(b, 2), "?")
            info["aid"] = _u16(b, 4) & 0x3FFF
            ie = b[6:]
        elif st == "auth" and len(b) >= 6:
            alg = _u16(b, 0)
            info.update(auth_algorithm=AUTH_ALG.get(alg, alg), auth_seq=_u16(b, 2), status=_u16(b, 4),
                        status_text=STATUS.get(_u16(b, 4), "?"))
            if alg == 3 and len(b) >= 8 and _u16(b, 2) == 1:
                info["sae_group"] = _u16(b, 6)
            ie = b[6:] if alg != 3 else None
        elif st in ("deauth", "disassoc") and len(b) >= 2:
            info["reason"], info["reason_text"] = _u16(b, 0), REASON.get(_u16(b, 0), "?")
        elif st in ("action", "action-no-ack") and len(b) >= 2:
            info["action"] = {"category": b[0], "category_name": ACTION_CAT.get(b[0] & 0x7F, b[0]), "code": b[1]}
            if b[0] == 127 and len(b) >= 8 and _mac(b[1:4]) == "18:fe:34":
                info["esp_now"] = _espnow(b, info)
            if b[0] == 4 and b[1] == 9 and len(b) >= 6 and _mac(b[2:5]) == "50:6f:9a":
                info["action"]["detail"] = "Wi-Fi Alliance public action (P2P/DPP/NAN)"
                rid = droneid.from_nan_action(b)                       # Remote ID over Wi-Fi NAN (research)
                if rid:
                    info["remote_id"] = rid
        if ie is not None:
            info.update(parse_elements(ie))
            if st != "probe-req" and info.get("ssid") == "<broadcast>":
                info["ssid"] = "<hidden>"                             # zero-length SSID in a beacon: hidden network
            sec = security(info)
            if sec:
                info["security"] = sec
            phy = phy_generation(info)
            if phy:
                info["phy"] = phy
    _vendors(info)
    _drone(info)
    info["summary"] = summary(info)
    return info


def _drone(info):
    """Flag the frame as a drone / drone controller, and say on what evidence.

    Remote ID is self-declared and authoritative (confidence 0.99 -- the aircraft said so). A drone
    vendor OUI or a drone-shaped SSID only identifies the platform, so those stay "-like" and below 1.
    """
    rid = info.get("remote_id")
    if rid:
        info["drone"] = {"platform": rid.get("model") or rid.get("uas_id") or "UAS",
                         "source": "Remote ID broadcast", "confidence": 0.99,
                         "evidence": droneid.summary(rid)}
        return
    w = info.get("wfb_ng")
    if w:
        info["drone"] = {"platform": "wfb-ng digital FPV / telemetry link (OpenIPC / OpenHD)",
                         "source": "wfb-ng 802.11 injection signature", "confidence": 0.9,
                         "evidence": f"link_id {w['link_id_hex']}" + (" (default link domain)" if w["default_link_domain"]
                                                                      else "") + f", radio port {w['radio_port']} ({w['stream']})"}
        return
    bp = info.get("elrs_backpack")
    if bp:
        info["drone"] = {"platform": f"ExpressLRS Backpack (bind UID {bp['bind_uid']})",
                         "source": "ESP-NOW MSP frame (CRC-valid)", "confidence": 0.85,
                         "evidence": f"MSP {bp['msp']}" + (f", VTX {bp['vtx_band']}{bp['vtx_channel']}"
                                                            if "vtx_band" in bp else "")}
        return
    mv = (info.get("ipv4") or {}).get("mavlink")
    if mv and mv.get("vehicle_type") is not None:
        info["drone"] = {"platform": f"MAVLink vehicle (sysid {mv.get('sysid')}, {mv.get('vehicle_type')})",
                         "source": "MAVLink telemetry over Wi-Fi (CRC-valid)", "confidence": 0.9,
                         "evidence": f"{mv['frames']} MAVLink frame(s): {', '.join(mv.get('messages', []))}"}
        return
    te = (info.get("ipv4") or {}).get("tello")
    if te:
        info["drone"] = {"platform": "DJI / Ryze Tello (Tello SDK)", "source": "Tello SDK over open Wi-Fi",
                         "confidence": 0.9 if te["port"] == "state" else 0.8,
                         "evidence": (f"state: battery {te.get('battery_pct')} %, height {te.get('height_cm')} cm"
                                      if te["port"] == "state" else f"command '{te['text']}' ({te['direction']})")}
        return
    tx = info.get("ta") or info.get("sa") or info.get("bssid")
    sig = droneid.signature(vendor(tx) if tx else None, info.get("ssid"))
    if sig:
        sig["source"] = "platform signature (no Remote ID in this frame)"
        sig["platform"] = sig["platform"] + "-like"
        info["drone"] = sig


def _caps(c):
    return [n for b, n in CAPS if c >> b & 1]


def _vendors(info):
    for k in ("addr1", "addr2", "addr3", "addr4", "addr1_ra", "addr2_ta", "sa", "da", "bssid", "ra", "ta"):
        if info.get(k):
            info[k + "_vendor"] = vendor(info[k])


def summary(info: dict) -> str:
    """One line for tables: what the frame is, who sent it to whom, and what it says."""
    st = info.get("subtype", "?")
    s = [st]
    if info.get("retry"):
        s.append("(retry)")
    if info.get("remote_id"):                        # a broadcasting aircraft outranks everything else
        return f"{st} " + droneid.summary(info["remote_id"])
    tx = info.get("sa") or info.get("ta")
    rx = info.get("da") or info.get("ra")
    if info.get("type") == "mgmt" and st in ("beacon", "probe-resp"):
        s.append(f"BSSID {info.get('bssid')} [{_short(info.get('bssid_vendor'))}]")
        s.append(f"SSID '{info.get('ssid', '?')}'")
        if info.get("ds_channel"):
            s.append(f"ch{info['ds_channel']}")
        for k in ("security", "phy"):
            if info.get(k):
                s.append(info[k])
        if info.get("ap_uptime"):
            s.append(f"up {info['ap_uptime']}")
        if info.get("bss_load"):
            s.append(f"{info['bss_load']['stations']} STA {info['bss_load']['channel_utilisation_pct']}% busy")
        w = info.get("wps") or {}
        dev = " ".join(x for x in (w.get("manufacturer"), w.get("model_name"), w.get("device_name")) if x)
        if dev:
            s.append(f"WPS '{dev}'")
    elif st == "probe-req":
        s.append(f"from {tx} [{_short(info.get('sa_vendor') or info.get('ta_vendor') or info.get('addr2_vendor'))}]")
        s.append(f"for '{info.get('ssid', '?')}'")
        if info.get("phy"):
            s.append(info["phy"])
        name = (info.get("p2p") or {}).get("device_name") or (info.get("wps") or {}).get("device_name")
        if name:
            s.append(f"device '{name}'")
    elif st in ("deauth", "disassoc"):
        s.append(f"{tx} -> {rx}: {info.get('reason_text')}")
    elif st == "auth":
        s.append(f"{tx} -> {rx}: {info.get('auth_algorithm')} seq {info.get('auth_seq')} {info.get('status_text')}")
    elif st.startswith(("assoc", "reassoc")):
        s.append(f"{tx} -> {rx}")
        if info.get("ssid"):
            s.append(f"SSID '{info['ssid']}'")
        if "status_text" in info:
            s.append(f"{info['status_text']} AID {info.get('aid')}")
    elif info.get("type") == "data":
        s.append(f"{info.get('direction', '')} {tx} -> {rx} BSSID {info.get('bssid')}")
        enc = info.get("encryption")
        if enc:
            s.append(f"{enc.get('cipher', 'encrypted')} PN {enc.get('pn')}" if enc.get("pn") is not None
                     else enc.get("cipher", "encrypted"))
        if info.get("eapol"):
            s.append(f"EAPOL {info['eapol'].get('handshake') or info['eapol'].get('eap', '')}")
        elif info.get("protocol"):
            s.append(info["protocol"])
        for k in ("ipv4", "ipv6"):
            if info.get(k):
                ip = info[k]
                s.append(f"{ip['src']} -> {ip['dst']} {ip.get('proto') or ip.get('next_header')}"
                         + (f" {ip.get('service')}" if ip.get("service") else "")
                         + (f" video {ip['video']['codec']} {ip['video']['nal']} ({ip['video']['transport']})"
                            if ip.get("video") else ""))
        if info.get("arp"):
            a = info["arp"]
            s.append(f"ARP {a['op']} {a['sender_ip']} -> {a['target_ip']}")
    elif info.get("type") == "ctrl":
        s.append(" -> ".join(x for x in (info.get("ta"), info.get("ra")) if x))
    elif st in ("action", "action-no-ack"):
        s.append(f"{tx} -> {rx}: {info['action']['category_name']} #{info['action']['code']}"
                 if info.get("action") else f"{tx} -> {rx}")
    if info.get("drone"):
        s.append(f"[{info['drone']['platform']}]")
    return " ".join(str(x) for x in s if x)


_SUFFIX = {"co", "co.", "ltd", "ltd.", "inc", "inc.", "corp", "corp.", "corporation", "pvt", "private", "limited",
           "llc", "gmbh", "company", "technology", "technologies", "systems", "communications", "electronics",
           "digital", "industrial", "international", "(india)", "india", "network", "networks", "co.,ltd",
           "co.,ltd.", "technologies,"}


def _short(v):
    if not v:
        return "?"
    if v.startswith("locally-administered"):
        return "randomized"
    words = v.split(",")[0].split()
    while len(words) > 1 and words[-1].lower() in _SUFFIX:
        words.pop()
    return " ".join(words[:3])
