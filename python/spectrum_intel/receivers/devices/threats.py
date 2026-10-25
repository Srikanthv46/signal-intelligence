"""Rules over the device tracks: what in this picture should an operator look at first.

Every finding states what was observed, what it means and how sure we are. Nothing here is a verdict:
a rule fires on evidence that is visible in the frames, and says so. Rules that cannot distinguish a
real attack from ordinary network behaviour say that too, at a lower severity.

    findings(devices, frames, links) -> list of findings, most serious first

Severity: critical (act now) | high | medium | low | info (inventory, not a problem).
Pure stdlib; no receiver imports.
"""
from __future__ import annotations

from collections import Counter, defaultdict

SEVERITY_ORDER = {"critical": 0, "high": 1, "medium": 2, "low": 3, "info": 4}

# A deauthentication burst is normal in small numbers (roaming, idle timeout). These are the points
# where the rate stops looking like housekeeping. Measured per capture, which is tens of ms, so even
# the lower bound is a high frame rate.
DEAUTH_SUSPECT = 8
DEAUTH_FLOOD = 25
WEAK_CIPHERS = ("WEP", "TKIP", "WPA1")
BEACON_FLOOD_APS = 15         # beaconing access points without an IEEE OUI ...
BEACON_FLOOD_PREFIXES = 8     # ... from this many different address prefixes


def _norm(mac):
    """MAC with the locally-administered bit cleared: virtual BSSIDs on one radio collapse onto the
    hardware address they were derived from."""
    return f"{int(mac[:2], 16) & ~0x02:02x}" + mac[2:]


def _same_radio(macs):
    """Do these BSSIDs look like virtual APs on one physical radio? Same OUI once the locally-
    administered bit is cleared, and at most one differing byte in the device part."""
    norm = [_norm(m) for m in macs]
    if len({n[:8] for n in norm}) != 1:
        return False
    tails = [n[9:].split(":") for n in norm]
    return all(sum(a != b for a, b in zip(tails[0], t)) <= 1 for t in tails[1:])


def _f(kind, severity, title, detail, *, confidence, evidence, devices=(), **extra):
    d = {"kind": kind, "severity": severity, "title": title, "detail": detail,
         "confidence": confidence, "evidence": evidence}
    if devices:
        d["devices"] = sorted(set(devices))
    d.update(extra)
    return d


# --------------------------------------------------------------------------------- rules
def _drones(devices):
    out = []
    for d in devices:
        dr = d.get("drone")
        if not dr:
            continue
        rid = d.get("remote_id") or {}
        where = ""
        if rid.get("lat") is not None:
            where = f" at {rid['lat']:.6f},{rid['lon']:.6f}"
            if rid.get("height_m") is not None:
                where += f", {rid['height_m']:.0f} m above takeoff"
        op = ""
        if rid.get("operator_lat") is not None:
            op = f" Operator at {rid['operator_lat']:.6f},{rid['operator_lon']:.6f}."
        if rid:
            out.append(_f("drone_remote_id", "high",
                          f"Drone broadcasting Remote ID: {dr['platform']}",
                          f"An unmanned aircraft is identifying itself{where}.{op} "
                          "Remote ID is transmitted in clear by regulation, so these are the "
                          "aircraft's own declared values.",
                          confidence=0.99, evidence=dr.get("evidence", ""), devices=[d["mac"]],
                          remote_id=rid, uas_id=rid.get("uas_id"),
                          operator_location=([rid["operator_lat"], rid["operator_lon"]]
                                             if rid.get("operator_lat") is not None else None)))
        else:
            out.append(_f("drone_signature", "medium",
                          f"Possible drone link: {dr['platform']}",
                          "A transmitter matches a drone platform by manufacturer and/or network name. "
                          "No Remote ID was broadcast, so this identifies the equipment, not a flight.",
                          confidence=dr.get("confidence", 0.7), evidence=dr.get("evidence", ""),
                          devices=[d["mac"]]))
    return out


def _deauth(devices, frames):
    out = []
    by_tx = Counter()
    targets = defaultdict(set)
    for fr in frames:
        m = fr["mac"] or {}
        if m.get("subtype") not in ("deauth", "disassoc"):
            continue
        tx = m.get("ta") or m.get("sa") or m.get("addr2")
        if tx:
            by_tx[tx] += 1
            rx = m.get("da") or m.get("ra") or m.get("addr1")
            if rx:
                targets[tx].add(rx)
    for tx, n in by_tx.items():
        if n < DEAUTH_SUSPECT:
            continue
        broadcast = any(t.lower() == "ff:ff:ff:ff:ff:ff" for t in targets[tx])
        sev = "critical" if n >= DEAUTH_FLOOD or broadcast else "high"
        out.append(_f("deauth_flood", sev,
                      f"Deauthentication burst from {tx}",
                      f"{n} deauthenticate/disassociate frames in this capture, aimed at "
                      f"{len(targets[tx])} address(es)"
                      + (", including the broadcast address, which disconnects every client at once"
                         if broadcast else "")
                      + ". This is how a client is forced off a network, either to deny it service or "
                        "to make it reconnect so the handshake can be captured. A busy AP shedding "
                        "idle clients can also produce these, so check whether one address is the "
                        "source of all of them.",
                      confidence=0.75 if sev == "high" else 0.85,
                      evidence=f"{n} deauth/disassoc frames from one transmitter in a single capture",
                      devices=[tx], n_frames=n, targets=sorted(targets[tx])[:16]))
    return out


def _evil_twin(devices):
    out = []
    by_ssid = defaultdict(list)
    for d in devices:
        if d.get("role") == "access point" and d.get("ssid") and not d["ssid"].startswith("<"):
            by_ssid[d["ssid"]].append(d)
    for ssid, aps in by_ssid.items():
        if len(aps) < 2:
            continue
        secs = {a.get("security") for a in aps}
        vendors = {(a.get("vendor") or "?").split(",")[0] for a in aps}
        open_now = [a for a in aps if a.get("security") in ("open", "WEP")]
        differs = len(secs) > 1
        if _same_radio([a["mac"] for a in aps]):
            out.append(_f("multi_bssid", "info",
                          f"'{ssid}' is advertised on {len(aps)} BSSIDs of one radio",
                          "The addresses differ only in the locally-administered bit and a byte, so "
                          "this is one access point running several virtual BSSs (guest network, "
                          "band steering), not a second device.",
                          confidence=0.85,
                          evidence=f"BSSIDs {', '.join(a['mac'] for a in aps)} share an OUI once the "
                                   "locally-administered bit is cleared",
                          devices=[a["mac"] for a in aps], ssid=ssid))
            continue
        if not differs and len(vendors) <= 1:
            continue                                  # same SSID, same vendor, same security: a normal mesh
        sev = "high" if (differs and open_now) else "medium"
        out.append(_f("duplicate_ssid", sev,
                      f"Network name '{ssid}' advertised by {len(aps)} different radios",
                      "The same network name is coming from more than one BSSID with "
                      + ("different security settings" if differs else "different hardware vendors")
                      + ". That is what an impersonating access point looks like; it is also what a "
                        "multi-AP home mesh or an extender looks like. The telling sign is a copy "
                        "offering weaker security than the original."
                      + (" One copy is open or WEP." if open_now else ""),
                      confidence=0.55 if sev == "medium" else 0.7,
                      evidence=f"BSSIDs {', '.join(a['mac'] for a in aps)}; "
                               f"security {sorted(str(s) for s in secs)}; vendors {sorted(vendors)}",
                      devices=[a["mac"] for a in aps], ssid=ssid,
                      security_per_bssid={a["mac"]: a.get("security") for a in aps}))
    return out


def _weak_security(devices):
    out = []
    for d in devices:
        if d.get("role") != "access point":
            continue
        sec = d.get("security")
        if sec is None:
            continue
        if sec == "open":
            out.append(_f("open_network", "low", f"Open network '{d.get('ssid')}'",
                          "No link-layer encryption: everything this network carries is readable "
                          "by anyone in range, and any device may associate.",
                          confidence=0.95, evidence=f"no RSN/WPA element, privacy bit clear ({d['mac']})",
                          devices=[d["mac"]], ssid=d.get("ssid")))
        elif any(w in sec for w in WEAK_CIPHERS):
            # TKIP offered *alongside* CCMP is a compatibility mode, not a broken network: clients
            # that can do CCMP will. Only a network with no strong cipher at all is actually weak.
            strong = "CCMP" in sec or "GCMP" in sec
            if strong:
                out.append(_f("legacy_cipher_offered", "low",
                              f"'{d.get('ssid')}' still accepts TKIP",
                              f"Advertises {sec}. Modern clients will use CCMP, but leaving TKIP "
                              "enabled keeps a weaker option available and blocks 802.11n rates.",
                              confidence=0.85, evidence=f"security = {sec} ({d['mac']})",
                              devices=[d["mac"]], ssid=d.get("ssid")))
            else:
                out.append(_f("weak_encryption", "medium", f"Obsolete encryption on '{d.get('ssid')}'",
                              f"Advertises {sec} with no modern cipher. WEP and TKIP are broken and "
                              "WPA1 is withdrawn; this traffic should not be treated as protected.",
                              confidence=0.9, evidence=f"security = {sec} ({d['mac']})",
                              devices=[d["mac"]], ssid=d.get("ssid")))
        w = d.get("wps") or {}
        if w and w.get("state") == "configured" and w.get("selected_registrar"):
            out.append(_f("wps_active", "medium", f"WPS registrar active on '{d.get('ssid')}'",
                          "The access point is in WPS enrolment mode, which accepts new devices with "
                          "a PIN or button press. It should not normally be left on.",
                          confidence=0.8, evidence=f"WPS selected registrar set ({d['mac']})",
                          devices=[d["mac"]], ssid=d.get("ssid")))
    return out


def _surveillance_and_identity(devices):
    """Devices that name themselves, and the networks a client has joined before."""
    out = []
    for d in devices:
        probed = d.get("probed_ssids") or []
        if len(probed) >= 1 and d.get("role") != "access point":
            out.append(_f("client_history", "info",
                          f"Client {d['mac']} is searching for known networks",
                          "A station lists networks it has joined before when it probes for them. "
                          "That history says where the device has been.",
                          confidence=0.9,
                          evidence=f"probe requests for {', '.join(repr(s) for s in probed[:8])}",
                          devices=[d["mac"]], probed_ssids=probed,
                          device_name=d.get("device_name"), vendor=d.get("vendor")))
        if d.get("device_name"):
            out.append(_f("named_device", "info",
                          f"Device identifies itself as '{d['device_name']}'",
                          "The device broadcasts its own name (WPS or Wi-Fi Direct), which usually "
                          "gives the model and often the owner.",
                          confidence=0.95, evidence=f"{d['mac']} [{d.get('vendor')}]",
                          devices=[d["mac"]], device_name=d["device_name"]))
    return out


def _eapol(frames):
    out = []
    n = Counter()
    who = defaultdict(set)
    for fr in frames:
        m = fr["mac"] or {}
        ev = m.get("eapol")
        if not ev:
            continue
        h = ev.get("handshake") or ev.get("eap") or "EAPOL"
        n[h] += 1
        for k in ("ta", "sa", "da", "ra"):
            if m.get(k):
                who[h].add(m[k])
    if n:
        out.append(_f("key_exchange", "info", "WPA key exchange observed",
                      "A client is joining or rekeying. The handshake itself does not expose the "
                      "password, but its presence marks exactly when a device associated.",
                      confidence=0.95,
                      evidence="; ".join(f"{k} x{v}" for k, v in n.most_common()),
                      devices=sorted({a for s in who.values() for a in s})[:8],
                      messages=dict(n)))
    return out


def _beacon_flood(devices):
    """Many access points with unregistered addresses from unrelated prefixes in one capture: forged beacons
    (mdk3 / mdk4-style flooding), which hide real networks and confuse clients (research, PLAN W6)."""
    fake = [d for d in devices if d.get("role") == "access point" and d.get("frame_types", {}).get("beacon")
            and (d.get("vendor") or "unknown OUI") in ("unknown OUI", "locally-administered (randomized/private)")]
    prefixes = {_norm(d["mac"])[:8] for d in fake}
    if len(fake) < BEACON_FLOOD_APS or len(prefixes) < BEACON_FLOOD_PREFIXES:
        return []
    return [_f("beacon_flood", "high", f"Beacon flood: {len(fake)} unregistered access points",
               "Access points with addresses no manufacturer registered, from unrelated address prefixes, all "
               "beaconing in one capture. That is what a beacon-flooding tool produces; it hides the real "
               "networks and confuses clients scanning for them.",
               confidence=0.75,
               evidence=f"{len(fake)} beaconing BSSIDs without an IEEE OUI, {len(prefixes)} different prefixes",
               devices=[d["mac"] for d in fake][:16], n_access_points=len(fake))]


def _drone_behaviour(links):
    """Signals whose Wi-Fi traffic has the shape of a drone video link without any identity (research, PLAN W1:
    receivers/wifi/field.py)."""
    out = []
    for dl in links:
        if dl.get("basis") != "behaviour":
            continue
        b = dl["bssid"]
        out.append(_f("drone_link_behaviour", "medium",
                      f"Possible drone video link: {b}" + (f" '{dl['ssid']}'" if dl.get("ssid") else ""),
                      "One client receives a steady, heavy downlink from this access point and sends back only "
                      "small frames: how a Wi-Fi drone streams video to the pilot's phone. Nothing identifies the "
                      "equipment, and an IP camera or a phone hotspot streaming video looks the same in one block.",
                      confidence=dl.get("confidence", 0.55),
                      evidence=", ".join(k for k, v in (dl.get("behaviour") or {}).items() if v),
                      devices=[b] + list(dl.get("clients") or []), bssid=b))
    return out


def _randomization(links):
    return [_f("mac_rotation", "low",
               f"{g['n_addresses']} addresses appear to be one device",
               "The addresses differ but the probe-request signature (and transmitter clock, where "
               "measured) match, which is what MAC randomisation on one handset looks like. Counting "
               "these as separate devices would overstate how many are present.",
               confidence=g["confidence"], evidence=g["evidence"], devices=g["macs"],
               cfo_hz_median=g.get("cfo_hz_median"))
            for g in links]


def _inventory(devices):
    aps = [d for d in devices if d.get("role") == "access point"]
    stas = [d for d in devices if d.get("role") != "access point"]
    if not devices:
        return []
    named = sorted({d["ssid"] for d in aps if d.get("ssid") and not d["ssid"].startswith("<")})
    return [_f("inventory", "info", f"{len(aps)} access points and {len(stas)} client devices heard",
               "Every transmitter decoded in this capture, with its vendor, signal level and activity.",
               confidence=1.0,
               evidence=(", ".join(f"'{s}'" for s in named[:8]) or "no named networks")
                        + (f" and {len(named) - 8} more" if len(named) > 8 else ""),
               n_access_points=len(aps), n_stations=len(stas))]


def findings(devices: list[dict], frames: list[dict], links: list[dict] | None = None,
             drone_links: list[dict] | None = None) -> list[dict]:
    out = []
    out += _drones(devices)
    out += _drone_behaviour(drone_links or [])
    out += _beacon_flood(devices)
    out += _deauth(devices, frames)
    out += _evil_twin(devices)
    out += _weak_security(devices)
    out += _eapol(frames)
    out += _randomization(links or [])
    out += _surveillance_and_identity(devices)
    out += _inventory(devices)
    out.sort(key=lambda f: (SEVERITY_ORDER.get(f["severity"], 9), -f.get("confidence", 0)))
    return out
