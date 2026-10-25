"""
catalog.py -- which drone video / data-link system a blind CP-OFDM link most resembles (research catalog).

Research, not identification: built from published specifications and third-party measurements
(research/drone_transmitters/02_digital_video.md); validated on synthetic links only. Every candidate lists what
matched and what conflicted, so the report never guesses silently.

Each system: Tu (us) or None (numerology not public), CP fraction(s), occupied bandwidth modes (MHz), bands (MHz),
channel centres or a raster, transmit antennas (1 / 2 / None = unknown or model-dependent), frame period (ms).
"""
from __future__ import annotations

RACEBAND = (5658, 5695, 5732, 5769, 5806, 5843, 5880, 5917)
DJI_FPV_CH = (5660, 5695, 5735, 5770, 5805, 5839, 5878, 5914)          # DJI O1 25 Mbps set; Walksnail uses it too
AR8030_EXTRA = (5758, 5788, 5828) + tuple(5340 + 40 * k for k in range(8))
ISM24, ISM58, U51 = (2400, 2484), (5650, 5950), (5150, 5350)
U58L = (5340, 5650)                                                   # HDZero L1-L8 low band
HDZERO_L = tuple(5362 + 37 * k for k in range(8))

SYSTEMS = [
    dict(name="DJI O4 (Mini 4 Pro, Air 3/3S, Avata 2, Neo/Flip, O4 Air Unit Pro)", tu=33.33, cp=(0.0703,),
         bw=(9.0, 18.0, 36.0, 54.0), bands=(ISM24, U51, ISM58), raster="dji_x5", antennas=2, frame_ms=2.5,
         src="M-lib + DJI specs (2T4R drones, O4 Air Unit Pro 2T2R)"),
    dict(name="DJI O4 Air Unit (Lite, 1T1R)", tu=33.33, cp=(0.0703,), bw=(9.0, 18.0, 36.0, 54.0), bands=(U51, ISM58),
         raster="dji_x5", antennas=1, frame_ms=2.5, src="DJI O4 Air Unit specs (1T1R)"),
    dict(name="DJI O3 / O3+ (Air 2S, Mini 3/3 Pro, Mavic 3, Avata, O3 Air Unit)", tu=66.67, cp=(0.0703,),
         bw=(9.0, 18.0, 36.0), bands=(ISM24, ISM58), raster="dji_x5", antennas=2, frame_ms=None,
         src="M-3p (Air 2S: 2048 @ 30.72 MSps, CP 144/160, SFBC), O3 Air Unit channels 5768.5 / 5804.5 / 5839.5"),
    dict(name="DJI OcuSync 2 (Mini 2, Mini SE 2, Air 2, Mavic 2)", tu=66.67, cp=(0.0703, 0.0781), bw=(9.0, 18.0),
         bands=(ISM24, ISM58), raster=None, antennas=None, frame_ms=None, src="M-lib (Mini 2)"),
    dict(name="DJI OcuSync 1 (Mavic Pro)", tu=66.67, cp=(0.0703,), bw=(1.3, 9.0, 18.0), bands=(ISM24,),
         raster="1mhz", antennas=None, frame_ms=None, src="DJI SDK (1.4 / 10 / 20 MHz), 1 MHz raster"),
    dict(name="Herelink (CubePilot)", tu=66.67, cp=(0.0703,), bw=(18.0, 9.0), bands=(ISM24,), raster=None,
         antennas=None, frame_ms=5.0, src="M-lib (LTE PSS/SSS; the herelink path names it when its PSS is found)"),
    dict(name="Autel SkyLink (EVO II / Lite+ / Nano / Max)", tu=66.67, cp=(0.0703,), bw=(9.0, 18.0),
         bands=(ISM24, ISM58, (900, 930)), raster=None, antennas=None, frame_ms=None,
         src="vendor C-UAS white paper (600 sc); RFUAV (Nano ~20 MHz)"),
    dict(name="Walksnail Avatar HD (dual-antenna air unit: Pro / V2 / Moonlight / GT)", tu=91.43, cp=(0.125,),
         bw=(18.9,), bands=(ISM58, U51), raster=DJI_FPV_CH + RACEBAND, antennas=2, frame_ms=7.163,
         src="waveshark #226 (802.16e 20 MHz, 1729 sc, STBC matrix A, frame 7.163 ms)"),
    dict(name="Walksnail Avatar HD (single-antenna air unit: Mini 1S / Mini 1S Lite / Nano)", tu=91.43, cp=(0.125,),
         bw=(18.9,), bands=(ISM58, U51), raster=DJI_FPV_CH + RACEBAND, antennas=1, frame_ms=7.163,
         src="waveshark #226; Caddx: Mini 1S / Lite / Nano V3 single antenna"),
    dict(name="Artosyn AR8030 link (Walksnail Ascent, BETAFPV ArtLynk, HGLRC Draco, Skyzone SteadyDigital, StartRC VT5)",
         tu=91.43, cp=(0.125,), bw=(4.7, 9.4, 18.9), bands=(ISM58,), raster=RACEBAND + AR8030_EXTRA, antennas=None,
         frame_ms=None, src="waveshark #228 / fpvwiki (5 / 10 / 20 MHz, 19 centres)"),
    dict(name="DJI Lightbridge 2 (Phantom 4 / Inspire 2 / M600)", tu=91.43, cp=(0.125,), bw=(9.2, 9.6),
         bands=(ISM24, ISM58), raster=None, antennas=None, frame_ms=None, src="M-lib (one Phantom 4 capture)"),
    dict(name="HDZero / Shark Byte (Divimath)", tu=37.926, cp=(0.03125,), bw=(23.44, 17.0), bands=(ISM58, U58L),
         raster=RACEBAND + (5760, 5800, 5705, 5740) + HDZERO_L, antennas=1, frame_ms=16.67,
         src="measured on an HDZero Freestyle V2 (27 MHz mode: 1024 FFT at 27 MSps, CP 32, 889 sc, 60 fps frame); "
             "17 MHz mode and Shark Byte not measured"),
    dict(name="DJI O1 (DJI FPV system, Air Unit, Caddx Vista, Polar, Runcam Link)", tu=None, cp=(), bw=(18.0, 36.0),
         bands=(ISM58,), raster=DJI_FPV_CH, antennas=None, frame_ms=None, src="DJI FCC filings (20 / 40 MHz)"),
    dict(name="Cambium PMP 450 (fixed-wireless infrastructure, not a drone)", tu=22.86, cp=(0.0625,), bw=(18.5,),
         bands=(ISM58, U51), raster=None, antennas=None, frame_ms=2.5, src="M-lib"),
]


def _in_bands(f, bands):
    return any(lo <= f <= hi for lo, hi in bands)


def _raster(f, r):
    if r == "dji_x5":                                        # DJI x.5 MHz centres
        return abs((f - 0.5) - round(f - 0.5)) <= 0.15
    if r == "1mhz":
        return abs(f - round(f)) <= 0.15
    return min(abs(f - c) for c in r) <= 1.0


def candidates(tu_us, cp_frac, bw_mhz, carrier_mhz, tx_antennas=None, frame_ms=None, top=3):
    """Ranked systems for a measured link. tx_antennas: 1 / 2 / None (from tx_streams: 2 only for STBC / SFBC;
    '1_stream_like' gives 1 but cannot exclude CDD on two antennas). Score = matches - conflicts (Tu mismatch
    excludes a system with a known Tu)."""
    out = []
    for s in SYSTEMS:
        m, c = [], []
        if s["tu"] is not None:
            if abs(tu_us - s["tu"]) > 0.015 * s["tu"] + 0.1:
                continue
            m.append(f"Tu {s['tu']} us")
        if s["cp"] and cp_frac is not None:
            (m if min(abs(cp_frac - x) for x in s["cp"]) <= 0.008 else c).append(f"CP {cp_frac:.4f} Tu")
        if bw_mhz is not None:
            near = min(s["bw"], key=lambda b: abs(bw_mhz - b))
            (m if abs(bw_mhz - near) <= 0.12 * near else c).append(f"BW {bw_mhz:.1f} MHz (mode {near})")
        if carrier_mhz is not None:
            if _in_bands(carrier_mhz, s["bands"]):
                m.append("band")
            else:
                c.append("band")
            if s["raster"] is not None:
                (m if _raster(carrier_mhz, s["raster"]) else c).append("channel raster")
        if tx_antennas is not None and s["antennas"] is not None:
            if tx_antennas == 2 and s["antennas"] == 1:
                c.append("2 TX antennas measured, system has 1")
            elif tx_antennas == 2 and s["antennas"] == 2:
                m.append("2 TX antennas")
        if frame_ms is not None and s["frame_ms"] is not None and abs(frame_ms - s["frame_ms"]) <= 0.05 * s["frame_ms"]:
            m.append(f"frame {s['frame_ms']} ms")                # bonus only: the blind frame estimate can be off
        if s["tu"] is None and len(m) < 2:                       # numerology unknown: needs BW + raster
            continue
        out.append(dict(system=s["name"], score=len(m) - len(c), matched=m, conflicts=c, source=s["src"]))
    out.sort(key=lambda d: (-d["score"], len(d["conflicts"])))
    return out[:top]


# ---------------------------------------------------------------- model level (research, PLAN M3 / M4; file 11 §1-2)
# catalog system -> (family, models sharing that system's RF; file 11 decision tables). A system whose models cannot be
# told apart on air keeps them all, so the model stays unresolved.
MODELS = {
    "DJI O4 (Mini 4 Pro, Air 3/3S, Avata 2, Neo/Flip, O4 Air Unit Pro)":
        ("DJI O4", ("Mini 4 Pro", "Air 3", "Air 3S", "Avata 2", "Neo", "Flip", "Matrice 4", "O4 Air Unit Pro")),
    "DJI O4 Air Unit (Lite, 1T1R)": ("DJI O4", ("O4 Air Unit (Lite)",)),
    "DJI O3 / O3+ (Air 2S, Mini 3/3 Pro, Mavic 3, Avata, O3 Air Unit)":
        ("DJI O3", ("Air 2S", "Mini 3 Pro", "Mavic 3 / 3 Classic / 3 Pro", "Avata", "Matrice 30", "O3 Air Unit")),
    "DJI OcuSync 2 (Mini 2, Mini SE 2, Air 2, Mavic 2)":
        ("DJI OcuSync 2", ("Mini 2 / Mini 4K / Mini 2 SE", "Mini 3", "Mavic Air 2", "Mavic 2 Pro / Zoom")),
    "DJI OcuSync 1 (Mavic Pro)": ("DJI OcuSync 1", ("Mavic Pro",)),
    "Herelink (CubePilot)": ("Herelink", ("Herelink (v1.1 / v2: identical on air)",)),
    "Autel SkyLink (EVO II / Lite+ / Nano / Max)": ("Autel SkyLink", ("EVO II", "EVO Lite+", "EVO Nano", "EVO Max")),
    "Walksnail Avatar HD (dual-antenna air unit: Pro / V2 / Moonlight / GT)":
        ("Walksnail Avatar HD", ("Avatar HD V1 / Pro", "Avatar HD V2 (dual antenna)", "Avatar GT", "Moonlight")),
    "Walksnail Avatar HD (single-antenna air unit: Mini 1S / Mini 1S Lite / Nano)":
        ("Walksnail Avatar HD", ("Mini 1S", "Mini 1S Lite", "Nano V3", "Avatar HD V2 (single antenna)", "Avatar GT2")),
    "Artosyn AR8030 link (Walksnail Ascent, BETAFPV ArtLynk, HGLRC Draco, Skyzone SteadyDigital, StartRC VT5)":
        ("Artosyn AR8030", ("Walksnail Ascent", "BETAFPV ArtLynk", "HGLRC Draco", "Skyzone SteadyDigital",
                            "StartRC VT5")),
    "DJI Lightbridge 2 (Phantom 4 / Inspire 2 / M600)": ("DJI Lightbridge 2", ("Phantom 4 Pro", "Inspire 2", "Matrice 600")),
    "HDZero / Shark Byte (Divimath)": ("HDZero / Shark Byte", ("HDZero VTX (any: Race / Freestyle / Whoop / ECO / AIO)",
                                                               "Shark Byte TX5S.1 / TX5M.1")),
    "DJI O1 (DJI FPV system, Air Unit, Caddx Vista, Polar, Runcam Link)":
        ("DJI O1 (DJI FPV system)", ("DJI FPV Air Unit", "DJI FPV drone", "Caddx Vista", "RunCam Link")),
}
# HDZero VTX models (research file 11 §2.3): one TX antenna each, one PHY; they differ only in power ceiling
HDZERO_VTX = ("Race V1", "Race V2", "Race V3", "Freestyle V1", "Freestyle V2", "Whoop", "Whoop Lite", "ECO", "AIO5",
              "AIO15")


def hdzero_model(bandwidth_mode="27 MHz (normal)"):
    """Model block for a link the hdzero path identified by its pilots: the family is certain, the model is not."""
    from spectrum_intel.core.drone_models import model_block
    ev = [f"HDZero (Divimath DM6300) pilot structure, {bandwidth_mode} mode, measured on an HDZero Freestyle V2",
          "every HDZero VTX transmits the same PHY; only the power ceiling differs (research file 11 section 2.3)",
          "Shark Byte (the Divimath predecessor) is not measured and not listed"]
    return model_block("HDZero", [f"HDZero {m} VTX" for m in HDZERO_VTX], ev, "inferred",
                       "every HDZero VTX is identical on air (power only)")


O3_AIR_UNIT_MHZ, O4_AIR_UNIT_MHZ = (5804.5, 5839.5), (5789.5, 5814.5)      # 20 / 10 MHz centres not shared (file 11)


def model_candidates(system_candidates, carrier_mhz=None, bw_mhz=None, tx_antennas=None):
    """Model block for a blind OFDM link from its ranked system candidates: the best-scoring conflict-free systems'
    models, narrowed by the DJI O3 / O4 channel centres, the HDZero 17 MHz mode and the 1-TX-only units."""
    from spectrum_intel.core.drone_models import model_block
    ok = [c for c in system_candidates or [] if not c["conflicts"] and c["system"] in MODELS]
    if not ok:
        return None
    best = max(c["score"] for c in ok)
    syst = [c["system"] for c in ok if c["score"] == best]
    ev = [f"system {s} (research catalog)" for s in syst]
    fams, models = [], []
    for s in syst:
        fam, ms = MODELS[s]
        fams.append(fam)
        models += ms
    if carrier_mhz is not None and any(f in ("DJI O3", "DJI O4") for f in fams):
        if min(abs(carrier_mhz - c) for c in O4_AIR_UNIT_MHZ) <= 0.3:
            keep = [f == "DJI O4" for f in fams]
            ev.append(f"channel {carrier_mhz:.1f} MHz is an O4 centre (not used by O3)")
        elif min(abs(carrier_mhz - c) for c in O3_AIR_UNIT_MHZ) <= 0.3:
            keep = [f == "DJI O3" for f in fams]
            ev.append(f"channel {carrier_mhz:.1f} MHz is an O3 centre (not used by O4)")
        else:
            keep = None
        if keep and any(keep):
            fams = [f for f, k in zip(fams, keep) if k]
            models = [m for s, k in zip(syst, keep) if k for m in MODELS[s][1]]
    if "HDZero / Shark Byte" in fams and bw_mhz is not None and abs(bw_mhz - 17.0) <= 2.0:
        models = [m if not m.startswith("HDZero") else "HDZero VTX with the RunCam Nano 90 camera (Whoop / Lite / AIO)"
                  for m in models]
        ev.append("17 MHz mode: only with the RunCam Nano 90 camera")
    if tx_antennas == 2:
        ev.append("2 transmit antennas (STBC / SFBC)")
    fam = " / ".join(dict.fromkeys(fams))
    return model_block(fam, models, ev, "inferred")


# uplink / control hop signatures (research, PLAN M6 / D3; file 11 §1 RFUAV Table 4 and file 02 §4: one published unit
# per model, so low confidence): (model, brand, hop width MHz, on-time ms, interval ms or None = not used); width within
# 15 %, on-time 20 %, interval 10 %. "O4-class" spans the measured O4 range (4.4-7 MHz, 0.41 ms, 2-2.5 ms)
HOP_MODELS = (("DJI FPV", "DJI", 5.0, 0.64, 4.0), ("DJI Avata 2", "DJI", 7.0, 0.41, 2.0),
              ("DJI Mavic 3 Pro", "DJI", 4.78, 1.7, 6.0), ("DJI Mini 3", "DJI", 3.5, 0.56, 5.96),
              ("DJI Mini 4 Pro", "DJI", 6.54, 0.404, 2.5),
              ("DJI OcuSync 2 aircraft (Mini 2-class; model unresolved)", "DJI", 1.1, 0.51, None),
              ("Herelink", "Herelink", 2.96, 0.52, 5.16), ("Skydroid H12", "Skydroid", 6.06, 0.25, 2.97),
              ("Skydroid T10", "Skydroid", 5.27, 0.27, 2.91), ("Siyi MK15", "Siyi", 3.32, 1.55, 4.99),
              ("Siyi MK32", "Siyi", 7.23, 0.47, 5.01), ("Siyi FT24", "Siyi", 4.59, 1.61, None))
HOP_RANGES = (("DJI O4 aircraft (model unresolved)", "DJI", (4.4, 7.0), 0.41, (2.0, 2.5)),)


def hop_model(hop_bw_mhz, hop_us, period_ms):
    """Model block for a hopping OFDM uplink from its hop width, on-time and interval, or None."""
    from spectrum_intel.core.drone_models import model_block
    if not hop_bw_mhz or not hop_us:
        return None
    on = hop_us / 1e3
    hits = [(m, b) for m, b, w, t, p in HOP_MODELS
            if abs(hop_bw_mhz - w) <= 0.15 * w and abs(on - t) <= 0.2 * t
            and (p is None or (period_ms and abs(period_ms - p) <= 0.1 * p))]
    hits += [(m, b) for m, b, (w0, w1), t, (p0, p1) in HOP_RANGES
             if 0.85 * w0 <= hop_bw_mhz <= 1.15 * w1 and abs(on - t) <= 0.2 * t
             and period_ms and 0.9 * p0 <= period_ms <= 1.1 * p1]
    if not hits:
        return None
    names = [m for m, _ in hits]
    ev = [f"hops {hop_bw_mhz} MHz x {hop_us} us every {period_ms} ms match the published timing of "
          f"{', '.join(names)} (one measured unit per model: low confidence)"]
    return model_block(" / ".join(dict.fromkeys(b for _, b in hits)) + " (hop timing)", names, ev, "inferred")
