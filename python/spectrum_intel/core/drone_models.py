"""
drone_models.py -- model-level identification tables and the common `model` block (research, PLAN M0-M2;
research/drone_transmitters/11_model_level_identification.md).

Facts re-typed with their sources (no code copied):
- DJI DroneID product-type byte -> model (file 11 section 0.1: RUB-SysSec DroneSecurity and Kismet tables, which agree
  on every shared code; Mini 2 = 63 confirmed by real decodes). Unknown codes stay unknown, never the nearest name.
- ANSI/CTA-2063-A serial-number manufacturer codes (ICAO-assigned, first 4 characters) -> maker (file 11 section 0.2).
  Add-on broadcast modules name the module maker, not the airframe.

`model_block` is the one shape every receiver reports at model level: the family, the model when exactly one fits,
the candidate list otherwise, the evidence, and how it was obtained ("decoded" identity vs "inferred" from RF
features). It never picks one of several candidates.
"""
from __future__ import annotations

DJI_PRODUCT_TYPES = {
    1: "Inspire 1", 2: "Phantom 3 series", 3: "Phantom 3 series", 4: "Phantom 3 Standard", 5: "Matrice 100",
    6: "ACEONE", 7: "WooKong-M", 8: "NAZA", 9: "A2", 10: "A3", 11: "Phantom 4", 12: "MG-1", 14: "Matrice 600",
    15: "Phantom 3 4K", 16: "Mavic Pro", 17: "Inspire 2", 18: "Phantom 4 Pro", 20: "N2", 21: "Spark",
    23: "Matrice 600 Pro", 24: "Mavic Air", 25: "Matrice 200", 26: "Phantom 4 series", 27: "Phantom 4 Advanced",
    28: "Matrice 210", 30: "Matrice 210 RTK", 31: "A3 AG", 32: "MG-2", 34: "MG-1A", 35: "Phantom 4 RTK",
    36: "Phantom 4 Pro V2.0", 38: "MG-1P", 40: "MG-1P RTK", 41: "Mavic 2", 44: "Matrice 200 V2 series",
    51: "Mavic 2 Enterprise", 53: "Mavic Mini", 58: "Mavic Air 2", 59: "Phantom 4 Multispectral", 60: "Matrice 300 RTK",
    61: "DJI FPV", 63: "Mini 2", 64: "Agras T10", 65: "Agras T30", 66: "Air 2S", 67: "Matrice 30", 68: "Mavic 3",
    69: "Mavic 2 Enterprise Advanced", 70: "Mini SE",
    # file 13 §3: DJI's AeroScope type map (codes 1-70 equal the table above); not yet confirmed by a real decode
    72: "Agras T40", 73: "Mini 3 Pro", 75: "Avata", 76: "Inspire 3", 77: "Mavic 3 Enterprise", 78: "FlyCart 30",
    82: "Agras T25", 83: "Agras T50", 84: "Mavic 3 Pro", 86: "Mavic 3 Classic", 87: "Mini 3", 88: "Mini 2 SE",
    89: "Matrice 350 RTK", 90: "Air 3", 91: "Matrice 3D / 3TD", 93: "Mini 4 Pro", 94: "Avata 2", 95: "Agras T60",
    96: "Agras T25P", 99: "Matrice 4T / 4E", 100: "Matrice 4D / 4TD", 102: "Flip", 103: "Matrice 400", 104: "Neo",
    105: "Air 3S", 108: "Mavic 4 Pro", 109: "Agras T100", 112: "Mini 5 Pro", 124: "Neo 2", 127: "Avata 360",
    128: "Agras T100S / T70S", 130: "Matrice 5TD", 132: "Agras T55", 133: "Lito X1", 140: "Lito 1",
    166: "Agras T200", 167: "FlyCart 200", 173: "Mini 4K",
}
DJI_CODES_CONFIRMED = 70          # codes <= 70: public decoder tables + real decodes (Mini 2 = 63); above: DJI map only

# CTA-2063-A manufacturer code -> (maker, kind); kind "module" = add-on Remote ID broadcaster, not the airframe
MFR_CODES = {
    "1581": ("DJI", "aircraft"), "1748": ("Autel Robotics", "aircraft"), "1588": ("Parrot", "aircraft"),
    "1668": ("Skydio", "aircraft"), "1596": ("Dronetag", "module"), "1863": ("XAG", "aircraft"),
    "2104": ("Talos", "aircraft"), "1865": ("Holy Stone", "aircraft"), "1910": ("Potensic", "aircraft"),
    "1587": ("AgEagle / senseFly", "aircraft"), "1792": ("uAvionix", "module"), "1817": ("Freefly", "aircraft"),
    "1849": ("Spektrum", "module"), "1930": ("SwellPro", "aircraft"), "1972": ("Lumenier", "aircraft"),
    "1973": ("FrSky", "module"), "2071": ("NewBeeDrone", "aircraft"), "1892": ("SJRC", "aircraft"),
    "1875": ("Bwine", "aircraft"), "1710": ("Atom UAV", "aircraft"), "1902": ("Elsight", "module"),
}
# file 13 §1: codes in the serials makers declared to the FAA (Remote ID DOCs, 2026-10-22): self-declarations, not the
# ICAO register; a code can name the RID module or serial issuer rather than the airframe maker
_BM = {"1906", "1908", "1925", "1961", "1872", "1774", "1787", "1796", "1798", "1809", "1891", "1893", "2003", "2009",
       "2055"}
for _code, _maker in {
        "1595": "Sunflower Labs", "1605": "Pyka", "1705": "Pyka", "1607": "Wingcopter", "1620": "Aero Systems West",
        "1633": "Wing Aviation", "1663": "Wing Aviation", "1634": "ACSL", "1641": "Applied Aeronautics",
        "1644": "Microdrones", "1656": "Quantum-Systems", "1698": "RigiTech", "1701": "Percepto", "1718": "Censys",
        "1720": "Leica Geosystems", "1728": "Sony", "1744": "Wingtra", "1759": "ABZ Innovation",
        "1767": "Teledyne FLIR", "1774": "unknown (declared make 'UAS')", "1785": "Watts Innovations",
        "1786": "Amazon Prime Air", "1787": "BlueMark", "1788": "Vision Aerial", "1790": "Hylio",
        "1791": "Inspired Flight", "1795": "Zephyr Systems", "1796": "Drone Defence", "1798": "Sentera",
        "1801": "Drone Volt", "1807": "Event38", "1808": "Draganfly", "1809": "Aerobits (idME)", "1810": "Matternet",
        "1812": "Ascent AeroSystems", "1813": "Harris Aerial", "1816": "Field / Explorer", "1836": "Coretronic",
        "1839": "Teal Drones", "1840": "Voliro", "1843": "blueflite", "1860": "XiL", "1862": "DeltaQuad",
        "1866": "SmartDrone", "1867": "Hoverfly", "1869": "Ruko", "1871": "Camflite", "1872": "Veeniix",
        "1873": "Zipline", "1877": "Lucid Bots", "1879": "Flightwave", "1884": "Workhorse Aero", "1886": "Flytrex",
        "1889": "WISPR", "1891": "Zing Drone Solutions", "1893": "Pierce Aerospace", "1894": "Drone-Clone Xperts",
        "1895": "EXO", "1898": "Fotokite", "1900": "Tucok", "1901": "Meteodrone", "1904": "ShenZhen C-Fly",
        "1906": "Aurelia", "1907": "Anzu Robotics", "1908": "MicroFlight", "1911": "Cider Copter", "1914": "BRINC",
        "1922": "TRIPLEFINE", "1924": "Arcsky", "1925": "DIY-D", "1926": "Skyfront", "1928": "XM2 Labs",
        "1933": "Triad Drones", "1937": "PrecisionVision", "1939": "SPECTA", "1943": "LISHI", "1947": "Swoop Aero",
        "1949": "SKYFISH", "1961": "Phoenix UAS", "1968": "Light Quadcopter UAV", "1969": "SkySpecs", "1971": "FIMI",
        "1984": "Nokia", "1990": "FACEGLE", "1995": "Blackbird", "2003": "Holyton", "2009": "Cinnamon Drones",
        "2021": "SHANGHENG", "2030": "Telegrid", "2046": "Aurora Flight Sciences", "2047": "Airobotics",
        "2049": "Aeroo", "2054": "Global Drone", "2055": "Gleesfun", "2061": "Aerostar Dynamics", "2068": "Black Swift",
        "2078": "CHCNAV", "2095": "Zero Zero Robotics", "2098": "XPECE", "2105": "Flyby Robotics",
        "2108": "Trace / Vesper UAS", "2110": "DMR Technologies", "2111": "TopXGun", "2123": "Dunlevy Consulting",
        "2124": "Avol", "2125": "Exedy Drones", "2126": "Caelone", "2129": "Heilongjiang Huida",
        "2132": "Angel Aerial Systems", "2137": "Antigravity", "2146": "Contixo", "2148": "Sustainable Skylines",
        "2157": "Ceres Air", "2160": "Rainmaker", "2161": "Tofect", "2164": "Terra Drone", "2165": "A2Z",
        "2183": "SubUAS", "2187": "Amor Fati", "2199": "Voltair", "2218": "Manna"}.items():
    MFR_CODES.setdefault(_code, (_maker, "module" if _code in _BM else "aircraft"))

# file 13 §2: serial prefix -> model(s) from the FAA DOC ranges ("list" prefixes, inferred from the first / last 1000
# listed serials, are not verified); (prefix, models, verified, DOC)
SERIAL_MODELS = [
    ("1581F0UY", ("Phantom 4 RTK",), True, "RID000000752"), ("1581F11V", ("Phantom 4 Pro V2.0",), True, "RID000000759"),
    ("1581F0M6", ("Mavic 2 Zoom",), True, "RID000000758"), ("1581F163", ("Mavic 2 Pro",), True, "RID000000757"),
    ("1581F4GC", ("Mavic 2 Enterprise Advanced",), True, "RID000000753"),
    ("1581F1WN", ("Mavic Air 2",), True, "RID000000754"), ("1581F3N3", ("Mavic Air 2",), True, "RID000000755"),
    ("1581F37Q", ("DJI FPV",), True, "RID000000756"), ("1581F3YT", ("Air 2S",), True, "RID000000007"),
    ("1581F4XF", ("Mini 3 Pro",), True, "RID000000001"), ("1581F6L5", ("Mini 3 Pro",), True, "RID000000001"),
    ("1581F5YH", ("Mini 3",), True, "RID000000062"), ("1581F6Z9", ("Mini 4 Pro",), True, "RID000000713"),
    ("1581F9DE", ("Mini 5 Pro",), True, "RID000002581"), ("1581F6N8", ("Air 3",), True, "RID000000714"),
    ("1581F895", ("Air 3S",), True, "RID000002046"), ("1581F4QW", ("Avata",), True, "RID000000008"),
    ("1581F6W8", ("Avata 2",), True, "RID000001376"), ("1581FA8J", ("Avata 360",), True, "RID000002651"),
    ("1581FBLK", ("Avata 360",), True, "RID000002650"), ("1581F45T", ("Mavic 3",), True, "RID000000009"),
    ("1581F4QZ", ("Mavic 3 Cine",), True, "RID000000010"), ("1581F67P", ("Mavic 3 Classic",), True, "RID000000041"),
    ("1581F67Q", ("Mavic 3 Pro",), True, "RID000000204"), ("1581F6MK", ("Mavic 3 Pro Cine",), True, "RID000000205"),
    ("1581F8LQ", ("Mavic 4 Pro",), True, "RID000002375"), ("1581F986", ("Mavic 4 Pro",), True, "RID000002376"),
    ("1581F4Z4", ("Inspire 3",), True, "RID000000192"), ("1581F5FH", ("Mavic 3 Enterprise M3E",), True, "RID000000019"),
    ("1581F5FJ", ("Mavic 3 Enterprise M3T",), True, "RID000000020"),
    ("1581F5FK", ("Mavic 3 Enterprise M3M",), True, "RID000000057"), ("1581FBEX", ("M3TA",), True, "RID000002584"),
    ("1581F6QA", ("Matrice 3D",), True, "RID000001166"), ("1581F6Q8", ("Matrice 3TD",), True, "RID000001167"),
    ("1581F4BN", ("Matrice 30",), True, "RID000000011"), ("1581F5BK", ("Matrice 30T",), True, "RID000000012"),
    ("1581F5BL", ("Matrice 30 Dock",), True, "RID000000273"), ("1581F5BM", ("Matrice 30T Dock",), True, "RID000000274"),
    ("1581F62H", ("Matrice 30T (EU / KR)",), True, "RID000001136"), ("1581F1ZN", ("Matrice 300 RTK",), True, "RID000000021"),
    ("1581F6GK", ("Matrice 350 RTK",), True, "RID000000256"), ("1581F7FV", ("Matrice 4E",), True, "RID000002198"),
    ("1581F7K3", ("Matrice 4T",), True, "RID000002199"), ("1581F8HH", ("Matrice 4D",), True, "RID000002246"),
    ("1581F8HG", ("Matrice 4TD",), True, "RID000002247"), ("1581F8DB", ("Matrice 400",), True, "RID000002486"),
    ("1581F7C6", ("FlyCart 30",), True, "RID000000760"), ("1581FAN4", ("FlyCart 100",), True, "RID000002620"),
    ("1581F446", ("Agras T10", "Agras T30"), True, "RID000001163 / 038"),
    ("1581F574", ("Agras T20P", "Agras T40"), True, "RID000001181 / 039"),
    ("1581F6BU", ("Agras T50 / T25",), True, "RID000001594"), ("1581F836", ("Agras T25P",), True, "RID000002587"),
    ("1581F8ZX", ("Agras T70P",), True, "RID000002579"), ("1581F8ZL", ("Agras T100",), True, "RID000002576"),
    ("1748CH", ("Autel EVO II",), True, "RID000000030"), ("1748FEV2H", ("Autel EVO II V3",), True, "RID000001969"),
    ("1748CLTC", ("Autel EVO Lite", "Autel EVO Lite+"), True, "RID000000032"),
    ("1748FEC1LTC", ("Autel EVO Lite",), True, "RID000001086"), ("1748FEC2LTC", ("Autel EVO Lite",), True, "RID000001086"),
    ("1748FEV5LTE", ("Autel EVO Lite",), True, "RID000001086"), ("1748CNA", ("Autel EVO Nano",), False, "RID000000320"),
    ("1748FEC3NAC", ("Autel EVO Nano",), True, "RID000001038"), ("1748FEC4NAC", ("Autel EVO Nano",), True, "RID000001038"),
    ("1748FEV3HMA", ("Autel EVO Max",), True, "RID000000151"), ("1748FEV3HMM", ("Autel MDX-1",), True, "RID000002344"),
    ("1748FEV3HMP", ("Autel MDX-1",), True, "RID000002344"), ("1748FEV3HMQ", ("Autel MDX-1",), True, "RID000002344"),
    ("1748FEV4AL", ("Autel Alpha",), True, "RID000001332"), ("1748CJD1", ("Autel Dragonfish Standard",), True, "RID000000048"),
    ("1748CJD2", ("Autel Dragonfish Lite",), True, "RID000000046"), ("1748CJD3", ("Autel Dragonfish Pro",), True, "RID000000047"),
    ("1748FDF1JD1", ("Autel Dragonfish Standard",), True, "RID000000048"),
    ("1748FDF2JD2", ("Autel Dragonfish Lite",), True, "RID000000046"),
    ("1748FDF3JD3", ("Autel Dragonfish Pro",), True, "RID000000047"),
    ("1588E040452", ("Parrot ANAFI Ai",), True, "RID000000067"), ("1588E040505", ("Parrot ANAFI USA",), True, "RID000000068"),
    ("1588E040465", ("Parrot ANAFI UKR GOV",), True, "RID000002663"), ("1668BE10JA00", ("Skydio X2E",), True, "RID000000097"),
    ("1668BS20", ("Skydio 2+",), True, "RID000000099"), ("1668BR40", ("Skydio X10 family",), False, "RID000000749"),
    ("1668BR70", ("Skydio R10",), False, "RID000002670"), ("1907F7USE23B", ("Anzu Raptor",), False, "RID000001035"),
    ("1907F7UEE23B", ("Anzu Raptor T",), False, "RID000001036"), ("1910F14", ("Potensic ATOM 2",), True, "RID000002237"),
    ("1910F916", ("Potensic RID-916 module",), True, "RID000001125"), ("1865AHS600", ("Holy Stone HS600",), False, "RID000000327"),
    ("1865CHS600D", ("Holy Stone HS600D",), False, "RID000001175"), ("1865CHS720G0", ("Holy Stone HS720G",), False, "RID000000411"),
    ("1865CHS920", ("Holy Stone HS920",), False, "RID000001150"), ("18699F11", ("Ruko F11",), False, "RID000000287"),
    ("1869AF11P", ("Ruko F11PRO",), False, "RID000000288"), ("1869BF11P2", ("Ruko F11PRO 2",), False, "RID000000683"),
    ("1869BF11GM", ("Ruko F11GIM",), False, "RID000000521"), ("1869CF11GM2", ("Ruko F11GIM2",), False, "RID000000251"),
    ("1869AU11P", ("Ruko U11PRO",), False, "RID000000286"), ("1869AU11S", ("Ruko U11S",), False, "RID000000285"),
    ("1892AF110", ("SJRC F11 4K PRO",), False, "RID000000414"), ("1892CF11S23", ("SJRC F11s 4K PRO",), False, "RID000000394"),
    ("1892CF22A23", ("SJRC F22 4K PRO",), False, "RID000000492"), ("1892CF22B23", ("SJRC F22s 4K PRO",), False, "RID000000379"),
    ("1892CF5SA00", ("SJRC F5s PRO",), False, "RID000000473"), ("1892BF7A238", ("SJRC F7 4K PRO",), False, "RID000000413"),
    ("1875BF7GB2", ("Bwine F7GB2",), False, "RID000000418"), ("1875CF7GB3", ("Bwine F7GB3",), True, "RID000002596"),
    ("1875CF7GM2", ("Bwine F7GIM2",), False, "RID000001711"), ("2095FH05", ("Zero Zero Robotics ZZ-H-1-003",), True, "RID000002593"),
    ("2095FH04", ("Zero Zero Robotics ZZ-H-1-004",), True, "RID000002594"), ("17287200", ("Sony Airpeak S1",), True, "RID000000017"),
    ("1767CSRH", ("Teledyne FLIR SkyRanger / SkyRaider",), False, "RID000000191"), ("1930FD02", ("SwellPro FD2",), True, "RID000001501"),
    ("1930FD03", ("SwellPro FD3",), True, "RID000001762"), ("1930FARD", ("SwellPro A-RID module",), True, "RID000001321"),
    ("1863FPEG8218", ("XAG P100",), True, "RID000000162"), ("1863FPEG8219", ("XAG P100 Pro",), True, "RID000000163"),
    ("1863FPEG8223", ("XAG P150",), True, "RID000002240"), ("1863F251", ("XAG P150 Max",), True, "RID000002602"),
    ("1863F242", ("XAG P60",), True, "RID000002588"),
]
_SERIAL_CHARS = set("0123456789ABCDEFGHJKLMNPQRSTUVWXYZ")              # CTA-2063-A excludes O and I


def dji_product(code):
    """DroneID product-type byte -> model name, or None (unknown codes are never mapped to the nearest name)."""
    try:
        return DJI_PRODUCT_TYPES.get(int(code))
    except (TypeError, ValueError):
        return None


def serial_model(serial):
    """A CTA-2063-A serial -> dict(models, verified, doc, prefix) from the FAA DOC serial prefixes (longest match), or
    None. The serial is matched with its length character (position 4) included, as the DOCs declare it."""
    if not isinstance(serial, str):
        return None
    hit = max((r for r in SERIAL_MODELS if serial.startswith(r[0])), key=lambda r: len(r[0]), default=None)
    return None if hit is None else dict(models=list(hit[1]), verified=hit[2], doc=hit[3], prefix=hit[0])


def dji_product_source(code):
    """Where the product-type name comes from (codes > 70: DJI's AeroScope map, no real decode yet)."""
    return ("public decoder tables (RUB-SysSec / Kismet), real decodes" if int(code) <= DJI_CODES_CONFIRMED else
            "DJI AeroScope type map (not yet confirmed by a real decode)")


def serial_manufacturer(serial):
    """A CTA-2063-A serial (4-char code, 1 length char, up to 15 chars) -> dict(code, maker, kind) or None."""
    if not isinstance(serial, str) or len(serial) < 6 or not set(serial) <= _SERIAL_CHARS:
        return None
    n = int(serial[4], 16) if serial[4] in "123456789ABCDEF" else 0
    if n == 0 or len(serial) != 5 + n:
        return None
    hit = MFR_CODES.get(serial[:4])
    return None if hit is None else dict(code=serial[:4], maker=hit[0], kind=hit[1])


def model_block(family, candidates=(), evidence=(), basis="inferred", unresolved_reason=None):
    """The common model-level report: model = the single candidate, else None with the candidate list."""
    cands = list(dict.fromkeys(candidates))
    model = cands[0] if len(cands) == 1 else None
    if model and basis != "decoded" and not model.endswith("-like"):
        model += "-like"
    out = dict(family=family, model=model, candidates=cands, basis=basis, evidence=list(evidence))
    if model is None:
        out["unresolved_reason"] = unresolved_reason or ("no model candidates" if not cands else
                                                         "several models share these on-air features")
    return out


def remote_id_model(rid):
    """Model block from a decoded Remote ID (core.remote_id fields): the manufacturer from the CTA-2063-A serial;
    per-model serial ranges are not tabulated yet (file 11: FAA DOC list), so the model stays unresolved."""
    mf, uas, ua = rid.get("manufacturer"), rid.get("uas_id"), rid.get("ua_type")
    ev = [f"Remote ID UA type: {ua}"] if ua is not None else []
    if mf and mf["kind"] == "aircraft":
        ev.insert(0, f"Remote ID serial {uas}: CTA-2063-A manufacturer code {mf['code']} = {mf['maker']}")
        sm = serial_model(uas)
        if sm:
            ev.insert(1, f"serial prefix {sm['prefix']} -> {' / '.join(sm['models'])} (FAA DOC {sm['doc']}"
                         + ("" if sm["verified"] else "; prefix inferred from the listed serials, unverified") + ")")
            return model_block(mf["maker"], sm["models"], ev, "decoded" if sm["verified"] else "inferred",
                               "the DOC prefix is shared by several models")
        return model_block(mf["maker"], (), ev, "decoded", "no FAA DOC serial prefix for this serial")
    if mf:
        ev.insert(0, f"Remote ID serial {uas}: add-on broadcast module by {mf['maker']} (code {mf['code']})")
        return model_block("unknown airframe (Remote ID module)", (), ev, "decoded",
                           "an add-on module names the module maker, not the airframe")
    if uas:
        ev.insert(0, f"Remote ID ID {uas}: no known manufacturer code")
    return model_block("unknown (Remote ID)", (), ev, "decoded", "no manufacturer code in the Remote ID")
