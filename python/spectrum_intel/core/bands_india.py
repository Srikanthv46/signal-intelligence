"""India frequency-allocation prior, 70 MHz - 6 GHz (70-300 MHz: research, signal gaps S14,
research/signal_gaps/06_vhf_below_300.md; the government land-mobile rows there are unconfirmed).

Sources: National Frequency Allocation Plan 2022 (DoT/WPC), IMT band plans
assigned in the 2010-2024 auctions, GSR 1048(E) / 2018-2021 delicensing
notifications. Operator names are a prior only — PSS/MIB carry no PLMN, so
reports never assert an operator from frequency alone.

Each entry: (f_lo_hz, f_hi_hz, key, description, expected signal classes)
"""
from __future__ import annotations

MHz = 1e6

BANDS = [
    (70 * MHz, 74.8 * MHz, "govt_mobile_70", "Fixed/mobile (government, defence)", ["narrowband_fm", "pmr"]),
    (74.8 * MHz, 75.2 * MHz, "marker_75", "Aeronautical marker beacons 75 MHz (ILS outer / middle / inner)",
     ["marker_beacon"]),
    (75.2 * MHz, 87.5 * MHz, "govt_mobile_75", "Fixed/mobile (government)", ["narrowband_fm", "pmr"]),
    (87.5 * MHz, 108 * MHz, "fm_bcast", "FM broadcasting (AIR / private FM), RDS", ["fm_broadcast", "rds"]),
    (108 * MHz, 111.975 * MHz, "vor_ils_loc", "VOR (even tenths) and ILS localizer (odd tenths)", ["vor", "ils"]),
    (111.975 * MHz, 117.975 * MHz, "vor", "VOR (and GBAS VDB)", ["vor"]),
    (117.975 * MHz, 137 * MHz, "air_vhf", "Aeronautical mobile: ATC AM voice 118-136.975, 121.5 emergency, ACARS, "
     "VDL Mode 2 136.725-136.975", ["am_voice", "acars", "vdl2"]),
    (137 * MHz, 138 * MHz, "metsat_137", "Meteorological satellite downlinks (Meteor-M LRPT)", ["metsat"]),
    (138 * MHz, 144 * MHz, "govt_mobile_138", "Fixed/mobile (defence)", ["narrowband_fm", "pmr"]),
    (144 * MHz, 146 * MHz, "ham_2m", "Amateur 2 m: FM repeaters, APRS, amateur satellites", ["aprs", "narrowband_fm"]),
    (146 * MHz, 156 * MHz, "govt_mobile_146", "Land mobile: police, state government, railways (FM, DMR)",
     ["narrowband_fm", "dmr", "pmr"]),
    (156 * MHz, 162.05 * MHz, "marine_vhf", "Maritime VHF: ch 16 156.8, DSC ch 70 156.525, AIS 161.975 / 162.025",
     ["marine_fm", "dsc", "ais"]),
    (162.05 * MHz, 174 * MHz, "govt_mobile_162", "Land mobile (police, railways, forest, PSU)",
     ["narrowband_fm", "dmr", "pmr"]),
    (174 * MHz, 225 * MHz, "vhf_band3", "Broadcasting Band III (VHF TV largely switched off; no DAB)", ["analog_tv"]),
    (225 * MHz, 300 * MHz, "mil_air_uhf", "Military aeronautical UHF (AM voice 225-400, 243.0 guard); fixed / mobile "
     "(defence)", ["am_voice"]),
    (300 * MHz, 328.6 * MHz, "govt_mobile_300", "Fixed/mobile (government, defence)", ["narrowband_fm", "pmr"]),
    (328.6 * MHz, 335.4 * MHz, "ils_glidepath", "ILS glide path (aeronautical)", ["ils"]),
    (335.4 * MHz, 380 * MHz, "govt_mobile_335", "Fixed/mobile (government, defence)", ["narrowband_fm", "pmr"]),
    (380 * MHz, 400 * MHz, "tetra_ppdr", "TETRA / PPDR public-safety trunking", ["tetra", "pmr"]),
    (400.05 * MHz, 406 * MHz, "metaids", "Meteorological aids (radiosondes) / Met-sat", ["radiosonde"]),
    (406 * MHz, 406.1 * MHz, "cospas_sarsat", "COSPAS-SARSAT distress beacons (EPIRB/PLB)", ["epirb"]),
    (406.1 * MHz, 433.05 * MHz, "govt_mobile_406", "Fixed/mobile (government), PMR", ["pmr", "narrowband_fm"]),
    (433.05 * MHz, 434.79 * MHz, "ism_433", "433 MHz SRD (delicensed, GSR 1048(E))", ["ism_srd", "lora"]),
    (434.79 * MHz, 450 * MHz, "govt_mobile_435", "Fixed/mobile (government), amateur 430-440", ["pmr", "narrowband_fm"]),
    (450 * MHz, 470 * MHz, "pmr_450", "PMR / land-mobile (police, utilities), DMR/analog FM", ["pmr", "dmr", "narrowband_fm"]),
    (470 * MHz, 585 * MHz, "uhf_tv", "UHF TV broadcasting (Doordarshan DVB-T2 / analog)", ["dvb_t2", "analog_tv"]),
    (585 * MHz, 698 * MHz, "uhf_bcast_imt", "Broadcasting / IMT (600 MHz, unassigned)", ["dvb_t2"]),
    (703 * MHz, 748 * MHz, "b28_ul", "LTE B28 / NR n28 uplink (FDD)", ["lte_ul", "nr_ul"]),
    (748 * MHz, 758 * MHz, "b28_gap", "B28 duplex gap", []),
    (758 * MHz, 803 * MHz, "b28_dl", "LTE B28 / NR n28 downlink (FDD) — Jio 5G SA, Airtel, Railways LTE-R", ["lte_dl", "nr_dl"]),
    (806 * MHz, 824 * MHz, "trunk_800", "Trunked radio / IMT 800 (B26 lower)", ["pmr", "lte_ul"]),
    (824 * MHz, 849 * MHz, "b5_ul", "LTE B5 / NR n5 uplink (850 MHz, FDD)", ["lte_ul"]),
    (851 * MHz, 869 * MHz, "trunk_850", "Trunked radio / B26 downlink lower", ["pmr", "lte_dl"]),
    (869 * MHz, 894 * MHz, "b5_dl", "LTE B5 / NR n5 downlink (850 MHz) — Jio, Airtel", ["lte_dl", "nr_dl"]),
    (865 * MHz, 867 * MHz, "ism_865", "865-867 MHz SRD (delicensed): LoRaWAN IN865, RFID", ["lora", "ism_srd", "rfid"]),
    (880 * MHz, 915 * MHz, "b8_ul", "GSM-900 / LTE B8 uplink (FDD)", ["gsm_ul", "lte_ul"]),
    (915 * MHz, 925 * MHz, "b8_gap", "900 MHz duplex gap / RFID (India 920-925? not delicensed)", []),
    (925 * MHz, 960 * MHz, "b8_dl", "GSM-900 / LTE B8 / NR n8 downlink — Airtel, Vi, BSNL", ["gsm_dl", "lte_dl", "nr_dl"]),
    (960 * MHz, 1164 * MHz, "arns_l", "Aeronautical radionavigation: DME/TACAN 962-1213, SSR 1030/1090, ADS-B", ["dme", "ssr", "adsb"]),
    (1164 * MHz, 1215 * MHz, "rnss_l5", "GNSS L5/E5a/E5b, NavIC L5 1176.45 (below noise) + DME", ["gnss", "dme"]),
    (1215 * MHz, 1300 * MHz, "rnss_l2_radar", "GNSS L2/G2/E6 + L-band surveillance radar", ["gnss", "radar"]),
    (1300 * MHz, 1350 * MHz, "arns_radar", "Aeronautical / L-band radar", ["radar"]),
    (1350 * MHz, 1400 * MHz, "fixed_radar_1350", "Fixed / mobile / radiolocation", ["radar"]),
    (1400 * MHz, 1427 * MHz, "passive_1400", "Passive (radio astronomy / EESS) — no emissions allowed", []),
    (1427 * MHz, 1518 * MHz, "l_band_imt", "L-band IMT / fixed links (not deployed for mobile)", ["fixed_link"]),
    (1518 * MHz, 1559 * MHz, "mss_dl", "MSS downlink (Inmarsat 1525-1559, Thuraya)", ["mss"]),
    (1559 * MHz, 1610 * MHz, "rnss_l1", "GNSS L1/E1/B1/G1 (below noise)", ["gnss"]),
    (1610 * MHz, 1626.5 * MHz, "mss_iridium", "MSS: Iridium 1616-1626.5 (TDMA bursts), Globalstar uplink", ["iridium"]),
    (1626.5 * MHz, 1660.5 * MHz, "mss_ul", "MSS uplink (Inmarsat terminals)", ["mss"]),
    (1660.5 * MHz, 1710 * MHz, "met_1675", "Met-sat / MetAids 1668-1710", ["metsat"]),
    (1710 * MHz, 1785 * MHz, "b3_ul", "LTE B3 / NR n3 uplink (1800 MHz FDD), GSM-1800 UL", ["lte_ul", "gsm_ul"]),
    (1785 * MHz, 1805 * MHz, "b3_gap", "1800 duplex gap", []),
    (1805 * MHz, 1880 * MHz, "b3_dl", "LTE B3 / NR n3 / GSM-1800 downlink — all operators", ["lte_dl", "nr_dl", "gsm_dl"]),
    (1880 * MHz, 1920 * MHz, "b39", "IMT 1900 (unassigned) / DECT-like cordless", []),
    (1920 * MHz, 1980 * MHz, "b1_ul", "LTE B1 / UMTS 2100 uplink", ["lte_ul", "umts_ul"]),
    (1980 * MHz, 2010 * MHz, "mss_s_ul", "MSS S-band uplink", []),
    (2010 * MHz, 2110 * MHz, "space_2025", "Space operation / research (uplink TT&C)", []),
    (2110 * MHz, 2170 * MHz, "b1_dl", "LTE B1 / NR n1 / UMTS 2100 downlink — Airtel, Vi, BSNL", ["lte_dl", "nr_dl", "umts_dl"]),
    (2170 * MHz, 2200 * MHz, "mss_s_dl", "MSS S-band downlink", []),
    (2200 * MHz, 2300 * MHz, "space_tlm", "Space research / telemetry (ISRO S-band TLM)", ["telemetry"]),
    (2300 * MHz, 2400 * MHz, "b40", "LTE B40 / NR n40 TDD (2300 MHz) — Jio, Airtel", ["lte_tdd", "nr_tdd"]),
    (2400 * MHz, 2483.5 * MHz, "ism_24", "2.4 GHz ISM (delicensed): Wi-Fi, BLE/BT, Zigbee, drones, microwave ovens",
     ["wifi", "wifi_b", "ble", "bt_classic", "zigbee", "dji", "fhss", "lora", "microwave_oven", "analog_video"]),
    (2483.5 * MHz, 2500 * MHz, "mss_24", "MSS downlink (Globalstar) / NavIC S-band 2492.028", []),
    (2500 * MHz, 2690 * MHz, "b41", "LTE B41 / NR n41 TDD (2500 MHz) — Vi, BSNL / BSS", ["lte_tdd", "nr_tdd"]),
    (2700 * MHz, 2900 * MHz, "radar_s_atc", "S-band ATC / weather radar", ["radar"]),
    (2900 * MHz, 3100 * MHz, "radar_s_marine", "Maritime / radiolocation radar", ["radar"]),
    (3100 * MHz, 3300 * MHz, "radar_s_def", "Radiolocation (defence radar)", ["radar"]),
    (3300 * MHz, 3670 * MHz, "n78", "NR n78 TDD (3.3-3.67 GHz) — Jio, Airtel, Vi 5G", ["nr_tdd"]),
    (3670 * MHz, 3700 * MHz, "c_band_low", "C-band (guard / FSS)", []),
    (3700 * MHz, 4200 * MHz, "fss_c", "C-band FSS satellite downlink (satellite TV; dish needed)", ["fss"]),
    (4200 * MHz, 4400 * MHz, "radalt", "Aircraft radio altimeters (FMCW)", ["radar_altimeter"]),
    (4400 * MHz, 4990 * MHz, "govt_4400", "Fixed/mobile (defence); 4.8-4.99 IMT-identified", ["fixed_link"]),
    (4990 * MHz, 5000 * MHz, "passive_4990", "Radio astronomy", []),
    (5000 * MHz, 5150 * MHz, "arns_mls", "ARNS / MLS 5030-5091 / FSS", []),
    (5150 * MHz, 5350 * MHz, "unii_12", "5 GHz Wi-Fi (U-NII-1/2, delicensed)", ["wifi", "dji", "analog_video"]),
    (5350 * MHz, 5470 * MHz, "radar_5400", "Radiolocation (no Wi-Fi)", ["radar"]),
    (5470 * MHz, 5725 * MHz, "unii_2e", "5 GHz Wi-Fi (U-NII-2e, delicensed) / weather radar 5600-5650", ["wifi", "radar"]),
    (5725 * MHz, 5875 * MHz, "ism_58", "5.8 GHz ISM / Wi-Fi U-NII-3 (delicensed): Wi-Fi, FPV video, DJI",
     ["wifi", "dji", "analog_video", "digital_vtx"]),
    (5875 * MHz, 5925 * MHz, "its_59", "ITS / 5.9 GHz", ["wifi"]),
    (5925 * MHz, 6425 * MHz, "wifi6e", "Lower 6 GHz (Wi-Fi 6E delicensing in progress in India)", ["wifi"]),
]


def lookup(f_hz: float) -> list[tuple]:
    """All allocations containing f (overlapping entries are allowed)."""
    return [b for b in BANDS if b[0] <= f_hz < b[1]]


def band_keys(f_lo: float, f_hi: float) -> list[str]:
    return [b[2] for b in BANDS if b[0] < f_hi and b[1] > f_lo]


def describe(f_hz: float) -> str:
    m = lookup(f_hz)
    return "; ".join(b[3] for b in m) if m else "no specific allocation in table"


# ---------------------------------------------------------------- channel rasters
def gsm_arfcn(f_dl_hz: float) -> int | None:
    f = f_dl_hz / MHz
    if 935 <= f <= 960:
        return int(round((f - 935) / 0.2))
    if 1805.2 <= f <= 1879.8:
        return int(round((f - 1805.2) / 0.2)) + 512
    return None


def lte_earfcn_dl(f_hz: float) -> tuple[int, int] | None:
    """(band, EARFCN) for an LTE downlink carrier centre (100 kHz raster)."""
    tab = {1: (2110, 0), 3: (1805, 1200), 5: (869, 2400), 8: (925, 3450), 28: (758, 9210),
           40: (2300, 38650), 41: (2496, 39650)}
    f = f_hz / MHz
    for band, (flo, noff) in tab.items():
        span = {1: 60, 3: 75, 5: 25, 8: 35, 28: 45, 40: 100, 41: 194}[band]
        if flo <= f <= flo + span:
            return band, int(round((f - flo) * 10)) + noff
    return None


def nr_arfcn(f_hz: float) -> int:
    """NR-ARFCN (global raster): 5 kHz steps < 3 GHz, 15 kHz steps 3-24.25 GHz."""
    f = f_hz / 1e3   # kHz
    if f < 3e6:
        return int(round(f / 5))
    return int(round(600000 + (f - 3e6) / 15))


WIFI_24 = {ch: 2407e6 + 5e6 * ch for ch in range(1, 14)}
WIFI_24[14] = 2484e6
WIFI_5 = {ch: 5000e6 + 5e6 * ch for ch in range(32, 178, 2)}
BLE_ADV = {37: 2402e6, 38: 2426e6, 39: 2480e6}
# allocations where radar / DME pulses are expected: classify r_pulsed; features.pulse_train runs only there
PULSED_BANDS = frozenset({"arns_l", "rnss_l5", "rnss_l2_radar", "arns_radar", "fixed_radar_1350", "radar_s_atc",
                          "radar_s_marine", "radar_s_def", "radar_5400", "unii_2e", "radalt"})


def pulsed_band(f_lo: float, f_hi: float, fc: float) -> bool:
    """The signal [f_lo, f_hi] (centre fc) overlaps a radar / DME allocation."""
    return bool((set(band_keys(f_lo, f_hi)) | {b[2] for b in lookup(fc)}) & PULSED_BANDS)


def nearest_wifi_channel(f_hz: float) -> tuple[int, float] | None:
    best = None
    for tab in (WIFI_24, WIFI_5):
        for ch, fc in tab.items():
            d = abs(f_hz - fc)
            if best is None or d < best[1]:
                best = (ch, d)
    return best


# frequency windows where FPV video links operate (1.2-1.3 GHz, 2.4 GHz, 5 GHz incl. the
# 4.9-5.3 GHz "lowrace"/custom range)
FPV_VIDEO_RANGES = [(1.06e9, 1.38e9), (2.33e9, 2.53e9), (4.90e9, 6.00e9)]
# analog VTX: the FPV ranges plus the 3.3 GHz analog band (3170-3470 MHz, two 8-channel groups; PLAN G9)
ANALOG_VTX_RANGES = [(1.06e9, 1.38e9), (2.33e9, 2.53e9), (3.15e9, 3.49e9), (4.90e9, 6.00e9)]
# blind digital (CP-OFDM) link search (PLAN G8): the FPV ranges widened to the professional / government UAV bands
# (L-band 1.35-1.444 GHz incl. the Chinese 1430-1444 MHz UAV band, 1.6-2.33 GHz, 4.4-4.9 GHz)
DIGITAL_LINK_RANGES = [(1.35e9, 1.444e9), (1.60e9, 2.53e9), (4.40e9, 6.00e9)]
_UAV_LINK_BANDS = [   # (lo MHz, hi MHz, context) for links found outside the ISM bands; first match wins
    (1430, 1444, "1.4 GHz UAV data-link band (China 1430-1444 MHz; Taisync-like)"),
    (1350, 1444, "L-band 1.35-1.44 GHz (government / professional UAV links)"),
    (1600, 2330, "1.6-2.33 GHz government / professional band (e.g. Skydio Connect MH, Microhard 1.8 GHz, "
                 "Silvus / Doodle Labs)"),
    (4400, 4900, "4.4-4.9 GHz government / professional band (e.g. Silvus, Doodle Labs, Persistent MPU5)"),
]


def uav_link_band(f_hz: float):
    """Band context of a link found outside the ISM / FPV bands (PLAN G8), else None."""
    f = f_hz / MHz
    return next((t for lo, hi, t in _UAV_LINK_BANDS if lo <= f < hi), None)


def cellular_at(f_hz: float) -> bool:
    """True inside an LTE / NR band of the table (keys b<n> / n<n>)."""
    return any(k[:1] in "bn" and k[1:2].isdigit() for k in band_keys(f_hz, f_hz))
