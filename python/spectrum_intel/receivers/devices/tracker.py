"""Turn decoded 802.11 frames into a per-device picture: who is there, what they are, what they do.

A single frame answers very little. The same address seen across a capture answers a lot: when it
appeared, how much it transmitted, who it talked to, how strong it was, which networks it has joined
before. That is what this builds.

    tracks(frames)      -> list of device dicts, strongest first
    fingerprint(mac)    -> what identifies a device when its MAC is randomised

Input is a list of normalised frame dicts (see `frame()`), so this module depends on no receiver: the
plugin feeds it whatever the 802.11b and OFDM receivers decoded, through the generic report contract.
Pure stdlib -- no numpy -- so it stays cheap next to the DSP.
"""
from __future__ import annotations

from collections import Counter

# Frame subtypes that only an access point sends.
AP_SUBTYPES = {"beacon", "probe-resp", "assoc-resp", "reassoc-resp"}
# Subtypes that carry no useful identity but do prove a transmitter exists.
CTRL_SUBTYPES = {"ack", "cts", "rts", "block-ack", "block-ack-req", "ps-poll", "cf-end"}
SEQ_GAP_MAX = 64          # a randomised address continuing another's 12-bit sequence counter within this many frames
MGMT_IDENTITY = {"beacon", "probe-resp", "probe-req", "assoc-req", "reassoc-req", "assoc-resp"}


def _is_group(mac):
    return bool(mac) and bool(int(mac[:2], 16) & 0x01)


def _is_random(mac):
    return bool(mac) and not int(mac[:2], 16) & 0x01 and bool(int(mac[:2], 16) & 0x02)


def frame(mac: dict, *, t0: float, dur_s: float = 0.0, snr_db=None, cfo_hz=None, channel=None,
          phy: str | None = None, rate_mbps=None, length_bytes=None, source: str = "") -> dict:
    """One decoded frame, normalised. `mac` is a mac80211.parse() result."""
    return {"mac": mac, "t0": t0, "dur_s": dur_s, "snr_db": snr_db, "cfo_hz": cfo_hz,
            "channel": channel, "phy": phy, "rate_mbps": rate_mbps,
            "length_bytes": length_bytes if length_bytes is not None else mac.get("length_bytes"),
            "source": source}


# --------------------------------------------------------------------------------- fingerprints
def _probe_fingerprint(m: dict) -> str | None:
    """What a probe request reveals about the chipset and OS, independent of the MAC address.

    The set and order of elements a station puts in a probe request is fixed by its driver, so it
    survives MAC randomisation. This is a grouping key, never an identity on its own.
    """
    if m.get("subtype") != "probe-req":
        return None
    ids = m.get("element_ids") or []
    if not ids:
        return None
    ext = m.get("ext_element_ids") or []
    vend = m.get("vendor_elements") or []
    parts = ["ie:" + ",".join(str(i) for i in ids)]
    if ext:
        parts.append("ext:" + ",".join(str(i) for i in ext))
    if vend:
        parts.append("vnd:" + ",".join(sorted(set(vend))))
    ht = m.get("ht") or {}
    if ht:
        parts.append(f"ht:{ht.get('spatial_streams')}ss{'/40' if ht.get('40mhz') else ''}")
    if m.get("rates_mbps"):
        parts.append("rates:" + ",".join(str(r) for r in m["rates_mbps"]))
    return "|".join(parts)


def _median(xs):
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def _round(v, n=1):
    return None if v is None else round(v, n)


# --------------------------------------------------------------------------------- the tracker
class _Device:
    def __init__(self, mac):
        self.mac = mac
        self.n_frames = 0
        self.t0 = self.t1 = None
        self.snr = []
        self.cfo = []
        self.subtypes = Counter()
        self.channels = Counter()
        self.phys = Counter()
        self.peers = Counter()
        self.probed = set()
        self.airtime = 0.0
        self.bytes = 0
        self.vendor = None
        self.is_ap = False
        self.ssid = None
        self.ssid_hidden = False
        self.bssid = None
        self.security = None
        self.phy = None
        self.uptime = None
        self.country = None
        self.bss_load = None
        self.wps = None
        self.device_name = None
        self.drone = None
        self.remote_id = None
        self.fingerprints = Counter()
        self.pn_max = None
        self.seqs = []
        self.seq_first = self.seq_last = None             # (time, sequence number) of the first / last frame
        self.rates = Counter()
        self.retries = 0
        self.pwr_save = 0
        self.sources = Counter()

    def see(self, f, m, role_ap):
        self.n_frames += 1
        t0 = f["t0"]
        self.t0 = t0 if self.t0 is None else min(self.t0, t0)
        self.t1 = max(self.t1 or t0, t0 + (f["dur_s"] or 0))
        self.airtime += f["dur_s"] or 0
        if f["length_bytes"]:
            self.bytes += f["length_bytes"]
        if f["snr_db"] is not None:
            self.snr.append(f["snr_db"])
        if f["cfo_hz"] is not None:
            self.cfo.append(f["cfo_hz"])
        if f["channel"] is not None:
            self.channels[f["channel"]] += 1
        self.subtypes[m.get("subtype", "?")] += 1
        if f["source"]:
            self.sources[f["source"]] += 1
        if f["rate_mbps"]:
            self.rates[f["rate_mbps"]] += 1
        if m.get("retry"):
            self.retries += 1
        if m.get("pwr_mgmt"):
            self.pwr_save += 1
        if m.get("seq") is not None:
            self.seqs.append(m["seq"])
            if self.seq_first is None or t0 < self.seq_first[0]:
                self.seq_first = (t0, m["seq"])
            if self.seq_last is None or t0 >= self.seq_last[0]:
                self.seq_last = (t0, m["seq"])
        self.is_ap = self.is_ap or role_ap
        self.vendor = self.vendor or m.get("ta_vendor") or m.get("sa_vendor") or m.get("bssid_vendor")
        fp = _probe_fingerprint(m)
        if fp:
            self.fingerprints[fp] += 1
        if m.get("subtype") == "probe-req" and m.get("ssid") not in (None, "<broadcast>"):
            self.probed.add(m["ssid"])
        name = (m.get("p2p") or {}).get("device_name") or (m.get("wps") or {}).get("device_name")
        if name:
            self.device_name = name
        if m.get("drone") and (self.drone is None or m["drone"]["confidence"] > self.drone["confidence"]):
            self.drone = m["drone"]
        if m.get("remote_id"):
            self.remote_id = m["remote_id"]
        if m.get("phy"):
            self.phys[m["phy"]] += 1
        enc = m.get("encryption") or {}
        if enc.get("pn") is not None:
            self.pn_max = max(self.pn_max or 0, enc["pn"])
        if role_ap:
            for k, a in (("ssid", "ssid"), ("security", "security"), ("phy", "phy"),
                         ("ap_uptime", "uptime"), ("country", "country"), ("bss_load", "bss_load"),
                         ("wps", "wps")):
                if m.get(k) is not None:
                    setattr(self, a, m[k])
            if m.get("ssid_hidden"):
                self.ssid_hidden = True

    def to_dict(self):
        d = {"mac": self.mac, "vendor": self.vendor,
             "address": "randomized (locally administered)" if _is_random(self.mac) else "global (OUI)",
             "role": "access point" if self.is_ap else "station",
             "n_frames": self.n_frames,
             "first_seen_s": _round(self.t0, 6), "last_seen_s": _round(self.t1, 6),
             "airtime_ms": _round(self.airtime * 1e3, 3),
             "bytes_seen": self.bytes or None,
             "snr_db_median": _round(_median(self.snr)), "snr_db_max": _round(max(self.snr) if self.snr else None),
             "cfo_hz_median": _round(_median(self.cfo), 0),
             "channels": dict(self.channels), "frame_types": dict(self.subtypes),
             "phy": self.phy or (self.phys.most_common(1)[0][0] if self.phys else None)}
        if self.is_ap:
            d.update(ssid="<hidden>" if (self.ssid_hidden and not self.ssid) else self.ssid,
                     ssid_hidden=self.ssid_hidden or None, security=self.security,
                     ap_uptime=self.uptime, country=(self.country or {}).get("code"),
                     bss_load=self.bss_load, wps=self.wps)
        else:
            d.update(associated_bssid=self.bssid, probed_ssids=sorted(self.probed) or None)
        if self.rates:
            d["rates_mbps"] = dict(self.rates)
        if self.retries:
            d["retry_fraction"] = round(self.retries / self.n_frames, 3)
        if self.pwr_save:
            d["power_save_fraction"] = round(self.pwr_save / self.n_frames, 3)
        if self.pn_max is not None:
            d["encryption_packet_counter_max"] = self.pn_max
        if self.device_name:
            d["device_name"] = self.device_name
        if self.peers:
            d["talks_to"] = [p for p, _ in self.peers.most_common(8)]
        if self.fingerprints:
            d["probe_fingerprint"] = self.fingerprints.most_common(1)[0][0]
        if self.seq_first is not None:
            d["seq_first"], d["seq_last"] = self.seq_first[1], self.seq_last[1]
        if self.drone:
            d["drone"] = self.drone
        if self.remote_id:
            d["remote_id"] = self.remote_id
        if self.sources:
            d["decoded_by"] = sorted(self.sources)
        return {k: v for k, v in d.items() if v is not None}


def tracks(frames: list[dict]) -> list[dict]:
    """Decoded frames -> one record per transmitting device, strongest first.

    A device is an address we have evidence transmitted: a TA/SA, or a BSSID that beacons. Addresses
    only ever seen as a destination are not devices we heard -- they are recorded as peers instead.
    """
    devs: dict[str, _Device] = {}

    def get(mac):
        if mac not in devs:
            devs[mac] = _Device(mac)
        return devs[mac]

    for f in frames:
        m = f["mac"] or {}
        st = m.get("subtype")
        tx = m.get("ta") or m.get("sa") or m.get("addr2")
        if st in AP_SUBTYPES and m.get("bssid") and not _is_group(m["bssid"]):
            tx = tx or m["bssid"]
        if not tx or _is_group(tx):
            continue
        d = get(tx)
        d.see(f, m, role_ap=(st in AP_SUBTYPES and m.get("bssid") == tx))
        rx = m.get("da") or m.get("ra") or m.get("addr1")
        if rx and not _is_group(rx):
            d.peers[rx] += 1
        b = m.get("bssid")
        if b and not _is_group(b) and b != tx and m.get("type") == "data":
            d.bssid = b

    # An address that beacons is an AP; anything that sends data to it is a station of that BSS.
    aps = {k for k, v in devs.items() if v.is_ap}
    for d in devs.values():
        if d.bssid is None and not d.is_ap:
            for p, _ in d.peers.most_common():
                if p in aps:
                    d.bssid = p
                    break

    out = [d.to_dict() for d in devs.values()]
    out.sort(key=lambda d: (-(d.get("snr_db_median") if d.get("snr_db_median") is not None else -99),
                            -d["n_frames"]))
    return out


def link_randomized(devices: list[dict], cfo_tol_hz: float = 2500.0) -> list[dict]:
    """Group randomised MACs that are probably one physical device.

    Two things survive a MAC change: the probe-request element signature (driver and OS) and the
    transmitter's crystal offset, which we already measure per frame. Agreement on both is evidence,
    not proof, so each group carries its confidence and the reason.
    """
    cand = [d for d in devices if d.get("address", "").startswith("randomized")
            and d.get("probe_fingerprint")]
    groups: list[dict] = []
    for d in cand:
        placed = False
        for g in groups:
            if g["probe_fingerprint"] != d["probe_fingerprint"]:
                continue
            if d.get("cfo_hz_median") is not None and g["cfo_hz_median"] is not None:
                if abs(d["cfo_hz_median"] - g["cfo_hz_median"]) > cfo_tol_hz:
                    continue
            g["macs"].append(d["mac"])
            g["n_frames"] += d["n_frames"]
            placed = True
            break
        if not placed:
            groups.append({"macs": [d["mac"]], "probe_fingerprint": d["probe_fingerprint"],
                           "cfo_hz_median": d.get("cfo_hz_median"), "n_frames": d["n_frames"]})
    by_mac = {d["mac"]: d for d in cand}
    out = []
    for g in groups:
        if len(g["macs"]) < 2:
            continue
        same_cfo = g["cfo_hz_median"] is not None
        # research, PLAN W7 (Vanhoef et al., AsiaCCS 2016): many drivers keep one sequence counter across a MAC
        # change, so the next address carrying on the previous one's counter ties the two together
        order = sorted(g["macs"], key=lambda a: (by_mac[a].get("first_seen_s") or 0.0, a))
        seq = any(by_mac[a].get("seq_last") is not None and by_mac[b].get("seq_first") is not None
                  and 0 < (by_mac[b]["seq_first"] - by_mac[a]["seq_last"]) % 4096 <= SEQ_GAP_MAX
                  for a, b in zip(order, order[1:]))
        out.append({"macs": sorted(g["macs"]), "n_addresses": len(g["macs"]), "n_frames": g["n_frames"],
                    "probe_fingerprint": g["probe_fingerprint"],
                    "cfo_hz_median": g["cfo_hz_median"],
                    "sequence_continues": seq,
                    "confidence": 0.9 if seq else 0.8 if same_cfo else 0.6,
                    "evidence": "identical probe-request element signature"
                                + (f" and transmitter clock offset within {cfo_tol_hz:.0f} Hz"
                                   if same_cfo else "")
                                + (" and one sequence-number counter continuing across the address change"
                                   if seq else "")
                                + " -- one device rotating its MAC address"})
    return sorted(out, key=lambda g: -g["n_frames"])
