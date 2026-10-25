"""Pipeline plugin: DJI OcuSync video / RC uplink / DroneID (wraps .detector)."""
from __future__ import annotations


from spectrum_intel.decode.emission_adapter import EmissionReceiver

DJI_BANDS = {"ism_24", "ism_58", "unii_12", "unii_2e"}


class DjiDecoder(EmissionReceiver):
    name = "dji"
    protocol = "DJI OcuSync/DroneID"
    key = "dji"
    bands = DJI_BANDS
    method = "receivers.dji (research receiver)"
    steps = ["15/30 kHz CP-OFDM fold on the 0.5 ms slot + DJI Zadoff-Chu root search (video downlink)",
             "RC uplink: 72-subcarrier ZC hops", "DroneID: 9-symbol 15.36 MSps burst, ZC sync, turbo decode + CRC"]

    def run_detector(self, ctx):
        from .detector import detect
        return detect(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc, {"decode_droneid": True})

    def pretty(self, top, em):
        like = all(e.features.get("evidence_tier") == "like" for e in em)   # inferred, not confirmed
        return {"droneid": "DJI DroneID", "uplink": "DJI OcuSync RC uplink"}.get(top, f"DJI OcuSync {top}") + \
            ("-like" if like else "")

    def summarise(self, em, sig):
        out = super().summarise(em, sig)
        st = [e.features.get("payload_state") for e in em if e.features.get("payload_state")]
        if st:                                       # research: CRC-24 / CRC-16 split (python-only decoder)
            out["fields"]["droneid_payload_state"] = "decoded" if "decoded" in st else st[0]
        fr = [e.features["droneid_frame"] for e in em if e.features.get("droneid_frame")]
        if fr:                                       # research (PLAN D4): 9-symbol OcuSync 2/3 vs 8-symbol legacy
            out["fields"]["droneid_frame"] = max(set(fr), key=lambda v: (fr.count(v), v))
        # model level (research, PLAN M1): the decoded product-type byte wins; the frame layout gives the generation
        from spectrum_intel.core.drone_models import dji_product, dji_product_source, model_block
        dt = next((e.features["device_type"] for e in em if e.features.get("device_type") is not None), None)
        if isinstance(dt, str) and not dt.strip().isdigit():
            # the external decoder (droneid_packet.DRONE_TYPES) already maps the byte to a name; "Unknown" = unmapped
            known = dt != "Unknown"
            out["fields"]["model"] = model_block(
                "DJI", [dt] if known else [], [f"DroneID product type: {dt} (the decoder's own table)"], "decoded",
                None if known else "product-type code not in the decoder's table")
        elif dt is not None:
            name = dji_product(dt)
            ev = [f"DroneID product type {dt}" + (f" = {name} ({dji_product_source(dt)})" if name
                                                  else " (code not in the public table)")]
            out["fields"]["model"] = model_block("DJI", [name] if name else [], ev, "decoded",
                                                 None if name else "product-type code not in the public table")
        elif fr:
            legacy = out["fields"]["droneid_frame"].startswith("8-symbol")
            # file 13 §4 (NDSS 2023, RUB legacy flag): Mavic Pro / Mavic 2 family; measured 2026-10-22 on 14 DJI Mini 4
            # Pro recordings (572 us bursts, first symbol absent, payload CRC-24 failing: O4 encryption)
            cands = ("Mavic Pro", "Mavic Pro Platinum", "Mavic 2 Pro", "Mavic 2 Zoom", "Mavic 2 Enterprise",
                     "Mavic 2 Enterprise Advanced", "Mini 4 Pro") if legacy else ()
            out["fields"]["model"] = model_block(
                "DJI", cands, [f"DroneID {out['fields']['droneid_frame']}"], "inferred",
                "the 8-symbol frame is sent by the Mavic Pro / Mavic 2 family and by O4 aircraft (Mini 4 Pro measured; "
                "other O4 models not measured)" if legacy else
                "the DroneID frame layout gives the generation only (payload not decoded)")
        return out

    def modulation(self, em):
        scs = em[0].features.get("scs_khz")
        return f"OFDM{f' {scs} kHz SCS' if scs else ''} (ZC sync)"


PLUGINS = [DjiDecoder()]
