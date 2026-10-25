"""
models.py -- model level for RC and telemetry links (research, PLAN M5;
research/drone_transmitters/11_model_level_identification.md sections 4-5). Plugin-side only: rules on what the
detector already decoded or measured; the result is the common model block (core.drone_models.model_block).

On air, most RC transmitters of one protocol differ only in power, so they stay unresolved; what can be told apart is
the hardware class (ELRS radio chip from the air mode), the firmware generation / region build (FrSky v1 / v2 / LBT,
Crossfire band plan, ELRS domain) and, for SiK, the RFD 3.x air-rate set.
"""
from __future__ import annotations

from spectrum_intel.core.drone_models import model_block

SIK_STD_KBPS = {2, 4, 8, 16, 19, 24, 32, 48, 64, 96, 128, 192, 250}
SIK_RFD_KBPS = {12, 56, 100, 125, 200, 224, 500, 750}            # RFD900x / ux / 868x (3.x firmware) [U inference]


def _decoded(em, prefix):
    return [e.decoded for e in em if e.label.startswith(prefix) and e.decoded]


def model(em):
    """Model block for one rc_link signal (its emissions), or None."""
    fr = _decoded(em, "frsky_")
    if fr:
        ok = [d for d in fr if d.get("crc_ok")] or fr
        var = sorted({d["variant"] for d in ok})
        ev = [f"FrSky {v}" + {"D16 v2": ": ACCST v2 firmware generation", "D16 LBT v1": ": EU-LBT region build",
                              "D16 v1": ": ACCST v1 FCC"}.get(v, "") for v in var]
        return model_block(f"FrSky {' / '.join(var)}", ("FrSky radio / internal module", "multiprotocol module"), ev,
                           "decoded", "no published on-air fingerprint separates FrSky radios and multiprotocol modules")
    ad = [e.decoded for e in em if e.label == "adsl" and e.decoded]
    if ad:
        cats = sorted({d["aircraft_category"] for d in ad if d.get("aircraft_category")})
        return model_block("ADS-L device", ("SkyEcho 2", "FLARM (ADS-L firmware)", "OGN tracker", "SoftRF"),
                           [f"ADS-L Traffic, category {', '.join(cats) or 'n/a'}, address table "
                            f"{', '.join(sorted({d['address_table'] for d in ad}))}"],
                           "inferred", "ADS-L packets do not name the device")
    fa = [e.decoded["fields"] for e in em if (e.decoded or {}).get("protocol") == "FANET" and e.decoded.get("crc_ok")
          and e.decoded.get("fields")]
    if fa:
        mf = sorted({d["manufacturer"] for d in fa})
        return model_block("FANET device", mf, [f"FANET manufacturer ID {', '.join(sorted({d['manufacturer_id'] for d in fa}))}"],
                           "decoded", "the manufacturer ID names the maker, not the device model")
    zb = [e.decoded for e in em if e.label.startswith("ieee802154") and e.decoded and e.decoded.get("crc_ok")]
    if zb:
        v = sorted({d["vendor"] for d in zb if d.get("vendor")})
        if v:
            return model_block("Digi XBee (802.15.4 / Zigbee)", ("XBee 802.15.4", "XBee Zigbee 3.0", "XBee DigiMesh 2.4"),
                               [f"64-bit address with the Digi OUI 00 13 A2 ({', '.join(v)})"], "decoded",
                               "the XBee 2.4 GHz variants share the PHY; the firmware is not visible in the MAC header")
        return model_block("IEEE 802.15.4 device", (), ["CRC-verified 802.15.4 frames, no vendor OUI in the addresses"],
                           "decoded", "short addresses carry no vendor")
    ty = [e.decoded for e in em if e.label == "toy_drone" and e.decoded]
    if ty:
        names = sorted({d["protocol"].rsplit(" (", 1)[0] for d in ty})
        return model_block("Toy drone (nRF24 / XN297)", names, [f"{d}" for d in sorted({d['protocol'] for d in ty})],
                           "decoded", "several toy protocols match these packets")
    xn = [e.decoded for e in em if e.label in ("bayang", "xn297_like") and e.decoded]
    if xn:
        ok = any(d.get("checksum_ok") for d in xn)
        return model_block("Bayang (XN297) toy drone" if ok else "XN297 toy drone", (),
                           ["XN297 prefix 0x710F55" + (", Bayang checksum verified" if ok else "")], "decoded",
                           "many toy-drone brands share the XN297 / Bayang protocol")
    rl = [e.decoded for e in em if e.label == "radiolink" and e.decoded]
    if rl:
        subs = sorted({d["sub_protocol"] for d in rl})
        cand = {"air": ("T8FB", "T8S"), "surface": ("RC4GS V3", "RC6GS V3", "RC8X"), "DumboRC": ("DumboRC transmitter",),
                "RC4G": ("RC4G",)}
        models = [m for s_ in subs for m in cand.get(s_, ())]
        return model_block("RadioLink", models, [f"RadioLink {s_} protocol (sync word / header)" for s_ in subs],
                           "inferred", "the transmitters of one RadioLink sub-protocol share the air protocol")
    cc = [e.decoded for e in em if e.label in ("hott", "hitec", "sfhss") and e.decoded]
    if cc:
        p = cc[0]["protocol"]
        if p.startswith("Graupner"):
            return model_block("Graupner HoTT", ("mz-12 / mz-16 / mz-24 / mz-32", "MC-series HoTT", "HoTT module"),
                               ["CC2500 MSK + FEC, sync D391, CRC-verified HoTT packet"], "decoded",
                               "every HoTT transmitter uses the same air protocol")
        if p.startswith("Hitec"):
            mm = ("Minima receiver link",) if any(d.get("variant") == "Minima" for d in cc) else \
                ("Aurora 9 / Flash 7 / Flash 8 (Optima)",)
            return model_block("Hitec", mm, ["CC2500 MSK, sync word = TX ID, CRC-verified"], "decoded",
                               "Optima transmitters share the air protocol")
        hints = sorted({d["model_hint"] for d in cc if d.get("model_hint")})
        return model_block("Futaba", hints or ("S-FHSS transmitter (T6J / T8J / T10J / T14SG / ...)",),
                           ["S-FHSS packet, model bytes " + ", ".join(sorted({d["model_bytes"] for d in cc}))],
                           "decoded" if hints else "inferred", "model bytes 3-4 not in the known list")
    hb = [e.decoded for e in em if e.label == "hubsan" and e.decoded]
    if hb:
        models = sorted({m for d in hb for m in d.get("models", [])})
        return model_block("Hubsan", models, [f"model constant bytes in {len(hb)} packet(s)"], "decoded",
                           "no data packet with a model's constant bytes")
    a1 = _decoded(em, "afhds1")
    if a1:
        return model_block("FlySky AFHDS 1", ("FS-CT6B", "FS-T6", "Turnigy 9X FlySky module", "multiprotocol module"),
                           ["A7105 ID code, 21-byte AFHDS 1 packet, hop table from the TX ID"], "decoded",
                           "every AFHDS 1 transmitter uses the same air protocol")
    af = _decoded(em, "afhds2a")
    if af:
        return model_block("FlySky AFHDS 2A", ("FS-i6", "FS-i6X", "other AFHDS 2A TX / multiprotocol module"),
                           ["A7105 ID code" + (", CRC-verified packets" if any(d.get("crc_ok") for d in af) else "")],
                           "decoded", "every AFHDS 2A transmitter uses the same air protocol")
    fl = [e.decoded for e in em if e.label == "flarm_like" and e.decoded]
    if fl:
        ok = [d for d in fl if d.get("radio_id")]
        ev = [f"FLARM radio ID {ok[0]['radio_id']} ({ok[0]['address_type']})"] if ok else \
            ["FLARM sync word; header not read reliably"]
        return model_block("FLARM device", (), ev, "decoded", "the device type is in the obfuscated payload (not decoded)")
    cf = [e.decoded for e in em if e.label == "crossfire_sync" and e.decoded]
    if cf:
        g = next((d["grid"] for d in cf if d.get("grid")), None)
        ev = ["Crossfire 150 Hz sync word"] + ([f"band plan {g['band']} (region build, not the model)"] if g else [])
        return model_block("TBS Crossfire", ("Micro TX V2", "Nano TX", "full-size TX", "Tango 2 / Mambo internal"), ev,
                           "inferred", "the same protocol on every Crossfire TX; only power differs")
    sk = [e.decoded for e in em if e.label == "sik_telemetry" and e.decoded]
    if sk:
        kb = round(sk[0].get("air_rate_bps", 0) / 1e3)
        if kb in SIK_RFD_KBPS:
            return model_block("SiK telemetry (RFD 3.x)", ("RFD900x", "RFD900ux", "RFD868x"),
                               [f"air rate {kb} kb/s is only in the RFD 3.x firmware set"], "inferred",
                               "RFD900x / ux / 868x differ only in the band edges")
        return model_block("SiK telemetry", ("100 mW SiK board (3DR / Holybro / mRo / HM-TRP)", "RFD900"),
                           [f"air rate {kb} kb/s (standard SiK set)"], "inferred",
                           "every SiK radio with this air rate looks the same on air")
    ml = [e for e in em if e.label.startswith("mlrs")]
    if ml:
        modes = sorted({m for e in ml for m in e.features.get("mlrs_modes") or []})
        chips = {"mlrs_19hz_7x": "SX1276 / LR1121 (SX127x-compatible SF6)", "mlrs_31hz_900": "SX126x / LLCC68 / LR1121",
                 "mlrs_19hz_900": "SX126x / LLCC68 / LR1121", "mlrs_50hz_2g4": "SX128x / LR1121",
                 "mlrs_31hz_2g4": "SX128x / LR1121", "mlrs_19hz_2g4": "SX128x / LR1121"}
        cands = list(dict.fromkeys(chips[m] + "-based mLRS module" for m in modes))
        ok = [e.decoded for e in ml if (e.decoded or {}).get("crc_ok")]
        ev = [f"mLRS air mode {', '.join(modes)}"] + ([f"FrameSyncWord {ok[0]['frame_sync_word']} (per binding)"]
                                                      if ok else [])
        return model_block("mLRS", cands, ev, "decoded" if ok else "inferred",
                           "the SX126x and SX127x SF6 framings differ by 2 symbols only")
    el = [e for e in em if "elrs" in e.label and e.features.get("modulation") == "lora"]
    if el:
        f = el[0].features
        sub = el[0].center_mhz < 1500.0
        ev = [f"LoRa SF{f.get('sf')} {f.get('bw_hz', 0) / 1e3:g} kHz {'sub-GHz' if sub else '2.4 GHz'}"]
        ch = f.get("elrs_channel")
        if ch:
            ev.append(f"domain {ch[0]} (regulatory build)")
        if sub and f.get("sf") == 5:
            ev.append("sub-GHz SF5 (200 Hz Full / 250 Hz) exists only on LR1121")
            return model_block("ExpressLRS", ("LR1121-based module",), ev, "inferred")
        cands = ("SX127x-based module", "LR1121-based module") if sub else ("SX1280-based module", "LR1121-based module")
        return model_block("ExpressLRS", cands, ev, "inferred", "this LoRa mode exists on both radio chips")
    eg = [e.decoded for e in em if e.label == "elrs_gfsk" and e.decoded]
    if eg:
        d = eg[0]
        return model_block("ExpressLRS", ("LR1121-based module",),
                           [f"300 kb/s GFSK SYNC packet verified (sync word {d['sync_word']} = UID5 / UID4)",
                            "the GFSK modes exist only on LR1121 hardware"], "decoded")
    labs = {e.label for e in em}
    if "elrs_flrc_like" in labs:
        return model_block("ExpressLRS-like (FLRC)", ("SX1280-based module",), ["650 kb/s FLRC-like rate"], "inferred")
    if "elrs_gfsk_like" in labs and all(e.center_mhz >= 2400.0 for e in em):
        return model_block("ExpressLRS-like (GFSK)", ("LR1121-based module",),
                           ["300 kb/s GFSK-like rate at 2.4 GHz (K1000 / DK modes: LR1121 only)"], "inferred")
    return None
