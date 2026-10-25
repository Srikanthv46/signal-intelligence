"""Lookup tables for Mode S / ADS-B fields: ICAO address blocks, operators, code meanings."""
from __future__ import annotations

# ICAO 24-bit address allocation by State (ICAO Annex 10 Vol III, Part I, Ch. 9). A subset: India,
# its neighbours, the Gulf/Asia carriers seen over India, and the large blocks elsewhere.
_BLOCKS = [
    (0x004000, 0x0043FF, "Zimbabwe"), (0x008000, 0x00FFFF, "South Africa"), (0x010000, 0x017FFF, "Egypt"),
    (0x018000, 0x01FFFF, "Libya"), (0x020000, 0x027FFF, "Morocco"), (0x040000, 0x040FFF, "Ethiopia"),
    (0x04C000, 0x04CFFF, "Kenya"), (0x05A000, 0x05A3FF, "Maldives"), (0x064000, 0x064FFF, "Nigeria"),
    (0x06A000, 0x06A3FF, "Qatar"), (0x0A0000, 0x0A7FFF, "Algeria"), (0x0D0000, 0x0D7FFF, "Mexico"),
    (0x100000, 0x1FFFFF, "Russia"), (0x300000, 0x33FFFF, "Italy"), (0x340000, 0x37FFFF, "Spain"),
    (0x380000, 0x3BFFFF, "France"), (0x3C0000, 0x3FFFFF, "Germany"), (0x400000, 0x43FFFF, "United Kingdom"),
    (0x440000, 0x447FFF, "Austria"), (0x448000, 0x44FFFF, "Belgium"), (0x450000, 0x457FFF, "Bulgaria"),
    (0x458000, 0x45FFFF, "Denmark"), (0x460000, 0x467FFF, "Finland"), (0x468000, 0x46FFFF, "Greece"),
    (0x470000, 0x477FFF, "Hungary"), (0x478000, 0x47FFFF, "Norway"), (0x480000, 0x487FFF, "Netherlands"),
    (0x488000, 0x48FFFF, "Poland"), (0x490000, 0x497FFF, "Portugal"), (0x498000, 0x49FFFF, "Czechia"),
    (0x4A0000, 0x4A7FFF, "Romania"), (0x4A8000, 0x4AFFFF, "Sweden"), (0x4B0000, 0x4B7FFF, "Switzerland"),
    (0x4B8000, 0x4BFFFF, "Turkey"), (0x4CA000, 0x4CAFFF, "Ireland"), (0x508000, 0x50FFFF, "Ukraine"),
    (0x680000, 0x6803FF, "Bhutan"), (0x700000, 0x700FFF, "Afghanistan"), (0x702000, 0x702FFF, "Bangladesh"),
    (0x704000, 0x704FFF, "Myanmar"), (0x706000, 0x706FFF, "Kuwait"), (0x708000, 0x708FFF, "Laos"),
    (0x70A000, 0x70AFFF, "Nepal"), (0x70C000, 0x70C3FF, "Oman"), (0x70E000, 0x70EFFF, "Cambodia"),
    (0x710000, 0x717FFF, "Saudi Arabia"), (0x718000, 0x71FFFF, "South Korea"), (0x720000, 0x727FFF, "North Korea"),
    (0x728000, 0x72FFFF, "Iraq"), (0x730000, 0x737FFF, "Iran"), (0x738000, 0x73FFFF, "Israel"),
    (0x740000, 0x747FFF, "Jordan"), (0x748000, 0x74FFFF, "Lebanon"), (0x750000, 0x757FFF, "Malaysia"),
    (0x758000, 0x75FFFF, "Philippines"), (0x760000, 0x767FFF, "Pakistan"), (0x768000, 0x76FFFF, "Singapore"),
    (0x770000, 0x777FFF, "Sri Lanka"), (0x778000, 0x77FFFF, "Syria"), (0x780000, 0x7BFFFF, "China"),
    (0x7C0000, 0x7FFFFF, "Australia"), (0x800000, 0x83FFFF, "India"), (0x840000, 0x87FFFF, "Japan"),
    (0x880000, 0x887FFF, "Thailand"), (0x888000, 0x88FFFF, "Viet Nam"), (0x890000, 0x890FFF, "Yemen"),
    (0x894000, 0x894FFF, "Bahrain"), (0x895000, 0x8953FF, "Brunei"), (0x896000, 0x896FFF, "United Arab Emirates"),
    (0x899000, 0x8997FF, "Taiwan"), (0x8A0000, 0x8A7FFF, "Indonesia"), (0xA00000, 0xAFFFFF, "United States"),
    (0xC00000, 0xC3FFFF, "Canada"), (0xC80000, 0xC87FFF, "New Zealand"), (0xE00000, 0xE3FFFF, "Argentina"),
    (0xE40000, 0xE7FFFF, "Brazil"), (0xE80000, 0xE80FFF, "Chile"),
]


def country(icao: str) -> str | None:
    a = int(icao, 16)
    for lo, hi, name in _BLOCKS:
        if lo <= a <= hi:
            return name
    return None


# ICAO operator designator (ICAO Doc 8585) -> (operator, IATA code). Indian carriers and the
# foreign carriers that fly to or over India.
OPERATORS = {
    "IGO": ("IndiGo", "6E"), "AIC": ("Air India", "AI"), "AXB": ("Air India Express", "IX"),
    "SEJ": ("SpiceJet", "SG"), "AKJ": ("Akasa Air", "QP"), "VTI": ("Vistara", "UK"),
    "IAD": ("AirAsia India", "I5"), "LLR": ("Alliance Air", "9I"), "FLG": ("Fly91", "IC"),
    "SDG": ("Star Air", "S5"), "BDA": ("Blue Dart Aviation", "BZ"), "QTR": ("Qatar Airways", "QR"),
    "UAE": ("Emirates", "EK"), "ETD": ("Etihad", "EY"), "FDB": ("flydubai", "FZ"), "ABY": ("Air Arabia", "G9"),
    "OMA": ("Oman Air", "WY"), "SVA": ("Saudia", "SV"), "GFA": ("Gulf Air", "GF"), "KAC": ("Kuwait Airways", "KU"),
    "SIA": ("Singapore Airlines", "SQ"), "MAS": ("Malaysia Airlines", "MH"), "THA": ("Thai Airways", "TG"),
    "ALK": ("SriLankan Airlines", "UL"), "BBC": ("Biman Bangladesh", "BG"), "RNA": ("Nepal Airlines", "RA"),
    "CPA": ("Cathay Pacific", "CX"), "BAW": ("British Airways", "BA"), "DLH": ("Lufthansa", "LH"),
    "AFR": ("Air France", "AF"), "KLM": ("KLM", "KL"), "THY": ("Turkish Airlines", "TK"),
    "UAL": ("United Airlines", "UA"), "AAL": ("American Airlines", "AA"), "FDX": ("FedEx", "FX"),
    "UPS": ("UPS", "5X"), "SWR": ("Swiss", "LX"), "VIR": ("Virgin Atlantic", "VS"),
    "AXM": ("AirAsia", "AK"), "TGW": ("Scoot", "TR"), "JAI": ("Jet Airways", "9W"),
}

CATEGORY = {  # (typecode, ME 'CA' field) -> emitter category, DO-260B 2.2.3.2.5.2
    (4, 0): "no category info", (4, 1): "light (< 15,500 lb)", (4, 2): "small (15,500-75,000 lb)",
    (4, 3): "large (75,000-300,000 lb)", (4, 4): "high-vortex large (B757)", (4, 5): "heavy (> 300,000 lb)",
    (4, 6): "high performance (> 5 g, > 400 kt)", (4, 7): "rotorcraft",
    (3, 0): "no category info", (3, 1): "glider / sailplane", (3, 2): "lighter-than-air",
    (3, 3): "parachutist / skydiver", (3, 4): "ultralight / hang-glider", (3, 5): "reserved",
    (3, 6): "unmanned aerial vehicle", (3, 7): "space / trans-atmospheric vehicle",
    (2, 0): "no category info", (2, 1): "surface emergency vehicle", (2, 2): "surface service vehicle",
    (2, 3): "point obstacle", (2, 4): "cluster obstacle", (2, 5): "line obstacle",
    (2, 6): "reserved", (2, 7): "reserved",
}

CAPABILITY = {  # CA field of DF11/DF17
    0: "level 1 transponder", 4: "level 2+ transponder, on ground", 5: "level 2+ transponder, airborne",
    6: "level 2+ transponder, ground or airborne", 7: "DR != 0 or FS = 2-5 (alert/SPI)",
}

CF = {  # DF18 control field
    0: "ADS-B ES/NT, ICAO address", 1: "ADS-B ES/NT, non-ICAO address", 2: "TIS-B fine format",
    3: "TIS-B coarse format", 4: "TIS-B / ADS-R management", 5: "TIS-B non-ICAO address",
    6: "ADS-R (rebroadcast of UAT)", 7: "reserved",
}

FLIGHT_STATUS = {
    0: "airborne", 1: "on ground", 2: "alert, airborne", 3: "alert, on ground",
    4: "alert + SPI", 5: "SPI (ident)", 6: "reserved", 7: "not assigned",
}

SURVEILLANCE_STATUS = {0: "no condition", 1: "permanent alert (emergency squawk)",
                       2: "temporary alert (squawk changed)", 3: "SPI (ident)"}

EMERGENCY = {0: "no emergency", 1: "general emergency", 2: "lifeguard / medical", 3: "minimum fuel",
             4: "no communications", 5: "unlawful interference", 6: "downed aircraft", 7: "reserved"}

ACAS_RI = {  # DF0/DF16 reply information
    0: "no operating ACAS", 2: "ACAS, resolution capability inhibited", 3: "ACAS, vertical-only RAs",
    4: "ACAS, vertical and horizontal RAs", 8: "no max airspeed data", 9: "max airspeed <= 75 kt",
    10: "max airspeed 75-150 kt", 11: "max airspeed 150-300 kt", 12: "max airspeed 300-600 kt",
    13: "max airspeed 600-1200 kt", 14: "max airspeed > 1200 kt",
}

ACAS_SL = {0: "ACAS inoperative", 1: "ACAS standby", 2: "TA only", 3: "TA/RA", 4: "TA/RA", 5: "TA/RA",
           6: "TA/RA", 7: "TA/RA"}

NACP_EPU_M = {11: 3, 10: 10, 9: 30, 8: 93, 7: 185, 6: 556, 5: 926, 4: 1852, 3: 3704, 2: 7408, 1: 18520}
