"""Pipeline plugin: device tracking and findings across every decoded 802.11 frame in a capture.

Runs after all receivers, reads their decoded frames off the generic report contract (it imports no
receiver), and adds report["devices"], report["device_groups"] and report["findings"], plus
signal.decoded.fields.device_macs so a signal can be traced back to the devices heard in it.
"""
from __future__ import annotations

from spectrum_intel.decode.base import Plugin

# Where decoded 802.11 frames surface today. Both Wi-Fi plugins publish "ppdus"; each PPDU carries
# "mac" (the dissected frame) and, for an A-MPDU, "macs" (every MPDU in it).
PPDU_KEYS = ("ppdus",)


def _frames_of(sig, name):
    """One signal's decoded PPDUs -> normalised frames. Handles the 802.11b and OFDM field names."""
    from . import tracker as T
    fields = ((sig.decoded or {}).get("fields") or {})
    out = []
    for key in PPDU_KEYS:
        for p in fields.get(key) or []:
            # Only FCS-verified frames become devices. Every receiver attaches "mac" solely when the
            # CRC-32 passed, and this makes that invariant explicit: a frame with bit errors would
            # otherwise invent an emitter that was never there.
            if p.get("fcs_ok") is False:
                continue
            macs = p.get("macs") or ([p["mac"]] if p.get("mac") else [])
            if not macs:
                continue
            if "start_s" in p:                                  # 802.11b DSSS records
                t0 = p["start_s"]
                dur = max(0.0, (p.get("end_s") or t0) - t0)
                cfo = p.get("cfo_hz")
                phy = "802.11b"
            else:                                               # OFDM records
                t0 = (p.get("start_us") or 0.0) / 1e6
                dur = (p.get("dur_us") or 0.0) / 1e6
                cfo = (p.get("cfo_khz") or 0.0) * 1e3
                phy = p.get("fmt")
            share = dur / len(macs) if macs else dur
            for m in macs:
                out.append(T.frame(m, t0=t0, dur_s=share, snr_db=p.get("snr_db"), cfo_hz=cfo,
                                   channel=p.get("channel"), phy=phy,
                                   rate_mbps=p.get("rate_mbps") or p.get("ht_rate_mbps"),
                                   source=name))
    return out


class DeviceTracking(Plugin):
    name = "devices"
    protocol = "device tracking and findings (802.11)"
    capture_level = True

    def capture_signals(self, ctx):
        return []                                   # defines no signals of its own

    def postprocess(self, ctx, signals):
        from . import threats as TH
        from . import tracker as T

        frames, per_signal = [], []
        for s in signals:
            dec = s.decoded or {}
            name = (dec.get("method") or "").split()[0] or "wifi"
            fr = _frames_of(s, name)
            per_signal.append((s, fr))
            frames.extend(fr)
        if not frames:
            return {}

        devices = T.tracks(frames)
        groups = T.link_randomized(devices)
        dlinks = [d for s in signals for d in (((s.decoded or {}).get("fields") or {}).get("drone_links") or [])]
        found = TH.findings(devices, frames, groups, dlinks)

        by_mac = {d["mac"] for d in devices}
        for s, fr in per_signal:
            macs = sorted({(f["mac"].get("ta") or f["mac"].get("sa") or f["mac"].get("addr2"))
                           for f in fr} & by_mac)
            if macs and isinstance(s.decoded, dict):
                s.decoded.setdefault("fields", {})["device_macs"] = macs

        out = {"devices": devices, "findings": found,
               "device_summary": {
                   "n_devices": len(devices),
                   "n_access_points": sum(1 for d in devices if d.get("role") == "access point"),
                   "n_stations": sum(1 for d in devices if d.get("role") != "access point"),
                   "n_randomized": sum(1 for d in devices if d.get("address", "").startswith("randomized")),
                   "n_frames": len(frames),
                   "method": "device tracking over decoded 802.11 frames (receivers.devices)",
                   "severity_counts": {k: sum(1 for f in found if f["severity"] == k)
                                       for k in ("critical", "high", "medium", "low", "info")
                                       if any(f["severity"] == k for f in found)}}}
        if groups:
            out["device_groups"] = groups
        return out

    def finalize_report(self, rep):
        """Point each finding and device at the signals it was heard in."""
        where = {}
        for s in rep.get("signals", []):
            for mac in (s.get("decoded", {}).get("fields", {}) or {}).get("device_macs", []):
                where.setdefault(mac, []).append(s["id"])
        for d in rep.get("devices", []):
            if d["mac"] in where:
                d["signal_ids"] = where[d["mac"]]
        for f in rep.get("findings", []):
            ids = sorted({i for m in f.get("devices", []) for i in where.get(m, [])})
            if ids:
                f["signal_ids"] = ids


PLUGINS = [DeviceTracking()]
