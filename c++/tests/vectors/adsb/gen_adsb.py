"""Reference vectors for receivers/adsb (C++ test test_rx_adsb.cpp).

pyModeS (GPL-3) is used here only as a black-box oracle, through the Python receiver: the C++
decoders are written from ICAO Annex 10 / Doc 9871 / DO-260B and checked against these outputs.

  demod.json + env_*.f32   modes.demod() on windows of the real_1090_commb envelope (|slice| at 8 MSps)
  fields.json              messages.fields() on random and structured frames of every DF / typecode /
                           Comm-B register (the frame dicts as demod() makes them, CRC/parity real)
  tracks.json              tracks.tracks() on synthetic aircraft: CPR even/odd + local (previous fix,
                           receiver), surface positions, implausible jumps, BDS50/60 resolution, met,
                           interrogation, flags; frames after tracks() (commb.register) included

Run: .venv/bin/python c++/tests/vectors/adsb/gen_adsb.py
"""
import copy
import json
import math
import random
from pathlib import Path

import numpy as np
import pyModeS as pms

from spectrum_intel.core.dsp import Channelizer
import sys

from spectrum_intel.receivers.adsb import messages, modes

trk = sys.modules["spectrum_intel.receivers.adsb.tracks"]

import os

REPO = Path(__file__).resolve().parents[4]
OUT = Path(os.environ.get("ADSB_OUT", Path(__file__).parent))
SCALE = int(os.environ.get("ADSB_SCALE", "1"))          # > 1: a bigger fuzz run (not committed)
R = random.Random(int(os.environ.get("ADSB_SEED", "20261010")))


# ------------------------------------------------------------------ bit helpers
def crc(bits):
    return modes._crc24(bits)


def to_bits(v, n):
    return [int(c) for c in format(v, f"0{n}b")]


def frame_from_bits(data_bits, df, addr_or_ic=0, t_s=0.0, snr=20.0):
    """Data bits (without parity) -> frame dict exactly as modes.demod() builds it."""
    p = crc(data_bits)
    if df not in (17, 18):
        p ^= addr_or_ic
    bits = data_bits + to_bits(p, 24)
    n = len(bits)
    res = modes._residual(bits)
    hexmsg = "%0*X" % (n // 4, int("".join(map(str, bits)), 2))
    m = {"t_s": t_s, "df": df, "bits": n, "hex": hexmsg, "parity_residual": f"{res:06X}", "snr_db": snr,
         "aa": "%06X" % int("".join(map(str, bits[8:32])), 2) if df in (11, 17, 18) else None}
    if df not in (11, 17, 18):
        m["confirmed"] = True
        m["aa"] = m["parity_residual"]
    else:
        m["confirmed"] = True
    return m


def setf(d, a, b, v):
    """1-based inclusive bit field [a, b] of the list d."""
    n = b - a + 1
    d[a - 1:b] = to_bits(v & ((1 << n) - 1), n)


def rbits(n):
    return [R.getrandbits(1) for _ in range(n)]


CHARS = "#ABCDEFGHIJKLMNOPQRSTUVWXYZ#####_###############0123456789######"


def char_code(c):
    return CHARS.index(c)


# ------------------------------------------------------------------ ME / MB builders
def es(df, tc, me_rest, ca=5, icao=None, t_s=0.0):
    """DF17/18 extended squitter: me_rest = ME bits 6..56 (51 bits)."""
    icao = R.getrandbits(24) if icao is None else icao
    d = to_bits(df, 5) + to_bits(ca, 3) + to_bits(icao, 24) + to_bits(tc, 5) + me_rest
    return frame_from_bits(d, df, t_s=t_s, snr=round(R.uniform(5, 40), 1))


def me_ident(cs=None):
    me = [0] * 51
    setf(me, 1, 3, R.randrange(8))
    if cs is None:
        valid = [i for i, c in enumerate(CHARS) if c != "#"]
        codes = [R.choice(valid) if R.random() < 0.9 else R.randrange(64) for _ in range(8)]
    else:
        codes = [char_code(c) for c in cs.ljust(8, "_")]
    for i, c in enumerate(codes):
        setf(me, 4 + 6 * i, 9 + 6 * i, c)
    return me


def alt12(ft):
    n = int(round((ft + 1000) / 25))
    b = to_bits(max(0, min(n, 2047)), 11)
    return b[:7] + [1] + b[7:]


def me_pos(alt_bits, odd, lat_cpr, lon_cpr, ss=None):
    me = [0] * 51
    setf(me, 1, 2, R.randrange(4) if ss is None else ss)       # ME 6-7 surveillance status
    me[2] = R.getrandbits(1)                                      # ME 8 NIC-B
    me[3:15] = alt_bits                                           # ME 9-20
    me[15] = R.getrandbits(1)                                     # ME 21 time
    me[16] = odd                                                  # ME 22 F
    setf(me, 18, 34, lat_cpr)
    setf(me, 35, 51, lon_cpr)
    return me


def me_surface(mov, trk_ok, trk7, odd, lat_cpr, lon_cpr):
    me = [0] * 51
    setf(me, 1, 7, mov)                                           # ME 6-12
    me[7] = trk_ok                                                # ME 13
    setf(me, 9, 15, trk7)                                         # ME 14-20
    me[15] = R.getrandbits(1)
    me[16] = odd
    setf(me, 18, 34, lat_cpr)
    setf(me, 35, 51, lon_cpr)
    return me


def me_vel(sub=None, gs=None, trk=None, vr=None):
    me = rbits(51)
    if sub is None:
        sub = R.choice([1, 1, 1, 2, 3, 3, 4, 0, 5])
    setf(me, 1, 3, sub)
    if gs is not None:
        vx, vy = gs * math.sin(math.radians(trk)), gs * math.cos(math.radians(trk))
        me[8] = int(vx < 0)
        setf(me, 10, 19, int(round(abs(vx))) + 1)
        me[19] = int(vy < 0)
        setf(me, 21, 30, int(round(abs(vy))) + 1)
        if vr is not None:
            me[32] = int(vr < 0)
            setf(me, 33, 41, int(round(abs(vr) / 64)) + 1)
    for a, b in ((10, 19), (21, 30), (33, 41), (45, 51)):         # zero / max fields sometimes
        u = R.random()
        if u < 0.08:
            setf(me, a, b, 0)
        elif u < 0.12:
            setf(me, a, b, (1 << (b - a + 1)) - 1)
        elif u < 0.16:
            setf(me, a, b, 1)
    return me


def me_random_tc(tc):
    if 1 <= tc <= 4:
        return me_ident()
    if 5 <= tc <= 8:
        mov = R.choice([0, 1, 2, 5, 8, 9, 12, 13, 20, 38, 39, 60, 93, 94, 100, 108, 109, 115, 123, 124, 125, 127,
                        R.randrange(128)])
        return me_surface(mov, R.getrandbits(1), R.randrange(128), R.getrandbits(1), R.getrandbits(17),
                          R.getrandbits(17))
    if 9 <= tc <= 18 or 20 <= tc <= 22:
        a = R.random()
        alt = [0] * 12 if a < 0.05 else rbits(12) if a < 0.5 else alt12(R.uniform(-1000, 50000))
        return me_pos(alt, R.getrandbits(1), R.getrandbits(17), R.getrandbits(17))
    if tc == 19:
        return me_vel()
    me = rbits(51)
    if tc == 28:
        setf(me, 1, 3, R.choice([0, 1, 1, 1, 2, 2, 3]))
        if R.random() < 0.3:
            setf(me, 4, 6, 0)
    elif tc == 29:
        setf(me, 1, 2, R.choice([1, 1, 1, 0, 2]))
    elif tc == 31:
        setf(me, 1, 3, R.choice([0, 0, 1, 1, 2]))
        setf(me, 36, 38, R.choice([0, 1, 2, 2, 2, 3]))
    return me


def ac13():
    u = R.random()
    if u < 0.05:
        return [0] * 13
    if u < 0.55:
        b = alt12(R.uniform(-1000, 45000))
        return b[:6] + [0] + b[6:]
    return rbits(13)


def surv(df, field13, rest_head, mb=None, addr=None):
    """DF0/4/5/16/20/21: rest_head = bits 6..19 (14 bits), field13 = bits 20..32, mb = 56 bits."""
    d = to_bits(df, 5) + rest_head + field13 + (mb or [])
    if df == 0 or df == 16:
        pass
    return frame_from_bits(d, df, R.getrandbits(24) if addr is None else addr, snr=round(R.uniform(5, 40), 1))


def df16_mv():
    if R.random() < 0.4:
        return to_bits(0x30, 8) + rbits(48)
    return rbits(56)


# --- Comm-B MB builders (Doc 9871 layouts); noisy on purpose so the inference sees every case
def noisy_status(d, sb, a, b, p_on=0.85, val=None):
    on = R.random() < p_on
    d[sb - 1] = int(on)
    if on:
        setf(d, a, b, val() if val else R.getrandbits(b - a + 1))
    elif R.random() < 0.1:
        setf(d, a, b, R.getrandbits(b - a + 1))
    else:
        setf(d, a, b, 0)


def maybe_zero(d, a, b, p=0.9):
    setf(d, a, b, 0 if R.random() < p else R.getrandbits(b - a + 1))


def mb40():
    d = [0] * 56
    noisy_status(d, 1, 2, 13, val=lambda: R.randrange(0, 2600))
    noisy_status(d, 14, 15, 26, val=lambda: R.randrange(0, 2600))
    noisy_status(d, 27, 28, 39, val=lambda: R.randrange(1800, 2400))
    maybe_zero(d, 40, 47)
    noisy_status(d, 48, 49, 51)
    maybe_zero(d, 52, 53)
    noisy_status(d, 54, 55, 56)
    return d


def twos(v, n):
    return v + (1 << n) if v < 0 else v


def mb50(roll=None, trk=None, gs=None, rate=None, tas=None):
    d = [0] * 56
    if roll is None:
        noisy_status(d, 1, 2, 11, val=lambda: twos(R.randrange(-120, 120), 10))
        noisy_status(d, 12, 13, 23, val=lambda: R.getrandbits(11))
        noisy_status(d, 24, 25, 34, val=lambda: R.randrange(0, 360))
        noisy_status(d, 35, 36, 45, val=lambda: twos(R.randrange(-80, 80), 10))
        noisy_status(d, 46, 47, 56, val=lambda: R.randrange(0, 360))
        return d
    d[0] = 1; setf(d, 2, 11, twos(int(round(roll * 256 / 45)), 10))
    t = trk if trk <= 180 else trk - 360
    d[11] = 1; setf(d, 13, 23, twos(int(round(t * 512 / 90)), 11))
    d[23] = 1; setf(d, 25, 34, int(round(gs / 2)))
    d[34] = 1; setf(d, 36, 45, twos(int(round(rate * 256 / 8)), 10))
    d[45] = 1; setf(d, 47, 56, int(round(tas / 2)))
    return d


def mb60(hdg=None, ias=None, mach=None, vrb=None, vri=None):
    d = [0] * 56
    if hdg is None:
        noisy_status(d, 1, 2, 12, val=lambda: R.getrandbits(11))
        noisy_status(d, 13, 14, 23, val=lambda: R.randrange(80, 400))
        noisy_status(d, 24, 25, 34, val=lambda: R.randrange(50, 230))
        noisy_status(d, 35, 36, 45, val=lambda: twos(R.randrange(-150, 150), 10))
        noisy_status(d, 46, 47, 56, val=lambda: twos(R.randrange(-150, 150), 10))
        return d
    h = hdg if hdg <= 180 else hdg - 360
    d[0] = 1; setf(d, 2, 12, twos(int(round(h * 512 / 90)), 11))
    if ias is not None:
        d[12] = 1; setf(d, 14, 23, int(round(ias)))
    if mach is not None:
        d[23] = 1; setf(d, 25, 34, int(round(mach / 0.004)))
    d[34] = 1; setf(d, 36, 45, twos(int(round(vrb / 32)), 10))
    d[45] = 1; setf(d, 47, 56, twos(int(round(vri / 32)), 10))
    return d


def mb44(wind=None, temp=None):
    d = [0] * 56
    setf(d, 1, 4, R.choice([0, 1, 2, 3, 4, 4, 5, 9]))
    if wind is None:
        noisy_status(d, 5, 6, 23, val=lambda: (R.randrange(0, 300) << 9) | R.getrandbits(9))
        setf(d, 24, 34, twos(R.randrange(-400, 300), 11) if R.random() < 0.9 else R.getrandbits(11))
    else:
        d[4] = 1; setf(d, 6, 14, wind[0]); setf(d, 15, 23, int(round(wind[1] * 256 / 180)) % 512)
        setf(d, 24, 34, twos(int(round(temp / 0.25)), 11))
    noisy_status(d, 35, 36, 46, p_on=0.4)
    noisy_status(d, 47, 48, 49, p_on=0.4)
    noisy_status(d, 50, 51, 56, p_on=0.4)
    return d


def mb45():
    d = [0] * 56
    for sb, a, b in ((1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12), (13, 14, 15)):
        noisy_status(d, sb, a, b, p_on=0.5)
    noisy_status(d, 16, 17, 26, val=lambda: twos(R.randrange(-320, 240), 10))
    noisy_status(d, 27, 28, 38, p_on=0.5)
    noisy_status(d, 39, 40, 51, p_on=0.5)
    maybe_zero(d, 52, 56)
    return d


def mb_head(byte, rest_len=48):
    return to_bits(byte, 8) + rbits(rest_len)


def mb10():
    d = to_bits(0x10, 8) + rbits(48)
    if R.random() < 0.7:
        setf(d, 10, 14, 0)
    return d


def mb17():
    d = rbits(24) + [0] * 32
    if R.random() < 0.7:
        d[6] = 1                                                   # BDS2,0 supported
    if R.random() < 0.3:
        setf(d, 25, 56, R.getrandbits(32) if R.random() < 0.5 else 1 << R.randrange(32))
    return d


def mb20(cs=None):
    d = to_bits(0x20, 8)
    if cs is None:
        valid = [i for i, c in enumerate(CHARS) if c != "#"]
        codes = [R.choice(valid) if R.random() < 0.93 else R.randrange(64) for _ in range(8)]
    else:
        codes = [char_code(c) for c in cs.ljust(8, "_")]
    for c in codes:
        d += to_bits(c, 6)
    return d


def mb30():
    d = to_bits(0x30, 8) + rbits(48)
    if R.random() < 0.6:
        setf(d, 29, 30, R.randrange(3))
    if R.random() < 0.6:
        setf(d, 16, 22, R.randrange(48))
    return d


def random_mb():
    k = R.random()
    if k < 0.03:
        return [0] * 56
    if k < 0.10:
        return rbits(56)
    return R.choice([mb10, mb17, mb20, mb30, mb40, mb40, mb44, mb45, mb50, mb50, mb60, mb60, mb60])()


AMBIG = []                                                       # MB fields inferred as BDS50,BDS60


def ambiguous_pool(n=60):
    tries = 0
    while len(AMBIG) < n and tries < 400000:
        tries += 1
        mb = mb50() if R.random() < 0.5 else mb60()
        f = surv(20 if R.random() < 0.5 else 21, ac13(), rbits(14), mb)
        if pms.bds.infer(f["hex"], mrar=True) == "BDS50,BDS60":
            AMBIG.append(mb)
    print(f"ambiguous BDS50/60 pool: {len(AMBIG)} after {tries} tries")


# ------------------------------------------------------------------ fields vectors
def random_frame():
    k = R.random()
    if k < 0.45:
        df = 17 if R.random() < 0.85 else 18
        tc = R.choice(list(range(0, 32)) + [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 18, 19, 19, 19, 19, 20, 22, 28, 28, 29,
                                             29, 29, 31, 31, 31])
        return es(df, tc, me_random_tc(tc), ca=R.randrange(8))
    if k < 0.52:
        d = to_bits(11, 5) + to_bits(R.randrange(8), 3) + rbits(24)
        ic = R.choice([0, 1, 5, 15, 16, 17, 31, 32, 63, 64, 100, 127, R.randrange(128)])
        return frame_from_bits(d, 11, ic, snr=round(R.uniform(5, 40), 1))
    df = R.choice([0, 4, 5, 16, 20, 20, 20, 21, 21, 21])
    head = rbits(14)
    field = ac13() if df in (0, 4, 16, 20) else rbits(13)
    mb = None
    if df == 16:
        mb = df16_mv()
    elif df in (20, 21):
        mb = random_mb()
    return surv(df, field, head, mb)


def gen_fields(n):
    frames = [random_frame() for _ in range(n)]
    for _ in range(6):                                             # empty MB, explicit corner cases
        frames.append(surv(20, [0] * 13, rbits(14), [0] * 56))
    frames.append(surv(21, rbits(13), rbits(14), mb20("IGO6123")))
    for mb in AMBIG[:20]:
        frames.append(surv(20, ac13(), rbits(14), mb))
    for i, f in enumerate(frames):
        f["t_s"] = R.randrange(0, 1 << 20) / modes.FS
    messages.fields(frames)
    return frames


# ------------------------------------------------------------------ tracks scenarios
NL = pms.common.cprNL


def cpr_encode(lat, lon, odd, surface=False):
    span = 90.0 if surface else 360.0
    dlat = span / (60 - odd)
    yz = math.floor(131072 * ((lat % dlat) / dlat) + 0.5)
    rlat = dlat * (yz / 131072 + math.floor(lat / dlat))
    nl = NL(rlat) - odd
    dlon = span / nl if nl > 0 else span
    xz = math.floor(131072 * ((lon % dlon) / dlon) + 0.5)
    return yz % 131072, xz % 131072


OPS = ["IGO", "AIC", "AXB", "UAE", "QTR", "SIA", "XYZ", "BAW"]


def scenario(k):
    t0 = 0.0 if R.random() < 0.4 else 1.79e9 + R.uniform(0, 1e6)
    n_ac = R.choice([1, 1, 2, 3])
    lat0 = R.choice([R.uniform(8, 32), R.uniform(-89.9, 89.9), R.uniform(86.4, 87.2), R.uniform(-0.01, 0.01),
                     R.uniform(-60, 60)])
    lon0 = R.choice([R.uniform(68, 97), R.uniform(-180, 180), R.uniform(179.5, 180), R.uniform(-180, -179.5)])
    ref = None
    u = R.random()
    if u < 0.35:
        ref = [round(lat0 + R.uniform(-0.8, 0.8), 6), round(lon0 + R.uniform(-0.8, 0.8), 6)]
    elif u < 0.6:
        ref = [round(lat0 + R.uniform(-0.8, 0.8), 6), round(lon0 + R.uniform(-0.8, 0.8), 6), R.choice([0, 12.5, 920])]
    elif u < 0.65:
        ref = [R.choice([12, 13]), R.choice([77, 78])]                # integer coordinates
    span = R.choice([0.03, 2.0, 12.0, 40.0, 90.0])
    use_ambig = R.random() < 0.3                                     # python raises KeyError on some
    frames = []
    for a in range(n_ac):
        icao = R.choice([0x800000 + R.getrandbits(18), R.getrandbits(24), 0x896000 + R.getrandbits(12)])
        lat, lon = lat0 + R.uniform(-0.5, 0.5), lon0 + R.uniform(-0.5, 0.5)
        lat = max(-89.95, min(89.95, lat))
        alt = R.choice([R.uniform(1000, 41000), R.uniform(-500, 2000), 37000.0])
        gs, track = R.uniform(120, 520), R.uniform(0, 360)
        tas = gs + R.uniform(-60, 60)
        hdg = (track + R.uniform(-10, 10)) % 360
        mach = tas / 600
        ias = tas * R.uniform(0.6, 0.95)
        vr = R.choice([0, 1500, -2000, R.uniform(-4000, 4000)])
        surface = R.random() < 0.15
        cs = R.choice(OPS) + str(R.randrange(1, 9999)) if R.random() < 0.8 else "N" + str(R.randrange(100, 999)) + "AB"
        squawk_emerg = R.random() < 0.1
        n = R.randrange(8, 50)
        ts = sorted(R.uniform(0, span) for _ in range(n))
        if R.random() < 0.2 and n > 3:
            ts[2] = ts[1]                                          # equal times
        odd = R.getrandbits(1)
        for t in ts:
            # aircraft moves along its track
            dlat = gs * t / 3600 / 60 * math.cos(math.radians(track))
            dlon = gs * t / 3600 / 60 * math.sin(math.radians(track)) / max(0.05, math.cos(math.radians(lat)))
            la, lo = max(-89.99, min(89.99, lat + dlat)), ((lon + dlon + 180) % 360) - 180
            if R.random() < 0.04:
                la, lo = la + R.uniform(-3, 3), lo + R.uniform(-3, 3)     # jump
                la = max(-89.99, min(89.99, la))
            kind = R.random()
            if surface and kind < 0.4:
                odd ^= 1
                yz, xz = cpr_encode(la, lo, odd, surface=True)
                f = es(17, R.choice([5, 6, 7, 8]), me_surface(R.randrange(1, 125), 1, R.randrange(128), odd, yz, xz),
                       icao=icao)
            elif kind < 0.40:
                odd = odd ^ 1 if R.random() < 0.85 else odd
                yz, xz = cpr_encode(la, lo, odd)
                if R.random() < 0.03:
                    yz, xz = R.getrandbits(17), R.getrandbits(17)
                tc = R.choice([11, 11, 11, 12, 13, 18, 20, 22])
                f = es(17, tc, me_pos(alt12(alt) if tc < 19 else to_bits(int(alt / 3.28084) & 0xFFF, 12), odd, yz, xz),
                       icao=icao)
            elif kind < 0.52:
                f = es(17, 19, me_vel(sub=1, gs=gs + R.uniform(-5, 5), trk=(track + R.uniform(-8, 8)) % 360, vr=vr),
                       icao=icao)
            elif kind < 0.57:
                f = es(17, R.choice([4, 4, 3, 1]), me_ident(cs if R.random() < 0.9 else "XYZ999"), icao=icao)
            elif kind < 0.62:
                sq = "7700" if squawk_emerg else R.choice(["2727", "1234", "2727"])
                idb = squawk_bits(sq)
                f = surv(5, idb, rbits(14), addr=icao)
            elif kind < 0.66:
                f = surv(4, ac_alt(alt), rbits(14), addr=icao)
            elif kind < 0.72:
                f = surv(20, ac_alt(alt), rbits(14), mb50(R.uniform(-30, 30), track, gs, R.uniform(-3, 3), tas), addr=icao)
            elif kind < 0.78:
                f = surv(R.choice([20, 21]), ac_alt(alt), rbits(14),
                         mb60(hdg, ias if R.random() < 0.8 else None, mach if R.random() < 0.8 else None, vr, vr),
                         addr=icao)
            elif kind < 0.81:
                f = surv(20, ac_alt(alt), rbits(14), R.choice(AMBIG) if use_ambig else mb50(), addr=icao)
            elif kind < 0.83:
                f = surv(20, ac_alt(alt), rbits(14), mb44((R.randrange(0, 120), R.uniform(0, 360)), R.uniform(-60, 20)),
                         addr=icao)
            elif kind < 0.85:
                f = surv(20, ac_alt(alt), rbits(14), mb20(cs if R.random() < 0.8 else "ABC123"), addr=icao)
            elif kind < 0.87:
                f = surv(R.choice([0, 16]), ac_alt(alt), rbits(14), df16_mv() if R.random() < 0.5 else None, addr=icao)
                if f["df"] == 16 and f["bits"] == 56:
                    continue
            elif kind < 0.89:
                d = to_bits(11, 5) + to_bits(5, 3) + to_bits(icao, 24)
                f = frame_from_bits(d, 11, R.choice([0, 3, 17, 40]))
            elif kind < 0.93:
                tc = R.choice([28, 29, 29, 31, 31])
                f = es(17, tc, me_random_tc(tc), icao=icao)
            else:
                f = random_frame()
                if R.random() < 0.7:
                    f = None
                if f is None:
                    continue
            f["t_s"] = t
            frames.append(f)
    frames = [f for f in frames if not (f["df"] in (0, 16) and f["bits"] == 112 and f["df"] == 0)]
    messages.fields(frames)
    for m in frames:
        m["t"] = t0 + m["t_s"]
    R.shuffle(frames)
    frames.sort(key=lambda m: m["t_s"]) if R.random() < 0.5 else None
    pyref = tuple(ref) if ref else None
    case = {"ref": ref}
    try:
        case["aircraft"] = trk.tracks(frames, pyref)
    except Exception as e:                                         # noqa: BLE001
        case["error"] = f"{type(e).__name__}: {e}"[:200]
    case["frames"] = frames
    return case


def squawk_bits(sq):
    a, b, c, d = (int(x) for x in sq)
    bit = lambda v, k: (v >> k) & 1                                  # noqa: E731
    # C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4
    return [bit(c, 0), bit(a, 0), bit(c, 1), bit(a, 1), bit(c, 2), bit(a, 2), 0, bit(b, 0), bit(d, 0), bit(b, 1),
            bit(d, 1), bit(b, 2), bit(d, 2)]


def ac_alt(ft):
    b = alt12(ft)
    return b[:6] + [0] + b[6:]


# ------------------------------------------------------------------ demod vectors
def gen_demod():
    meta = json.loads((REPO / "testdata/adsb/real_1090_commb.json").read_text())
    iq = np.fromfile(REPO / "testdata/adsb/real_1090_commb.cf32", np.complex64)
    ch = Channelizer(iq, meta["fs"], meta["fc"])
    y, _ = ch.slice(modes.F_ADSB, 2e6, modes.FS)
    env = np.abs(y)
    noise = np.median(env)
    thr = 6 * noise * modes.SPU / 2 * 2
    full = modes.demod(env, thr)
    env.astype(np.float32).tofile(OUT / "env_full.f32")
    wins = []
    s3 = int(round(full[2]["t_s"] * modes.FS))                      # DF21 + DF17 0.9 ms apart
    for s0, n in ((s3 - 1500, 12000), (len(env) - 3000, 3000), (0, 30000)):
        msgs = modes.demod(env[s0:s0 + n], thr)
        wins.append({"start": s0, "n": n, "msgs": msgs})
    return {"thr": float(thr), "n_env": len(env), "full": full, "windows": wins}


def gen_cpr(n=1500):
    """pyModeS CPR decoders (black box): local (airborne/surface) against a reference, global pairs,
    NL at and around the zone boundaries."""
    loc, glob, nl = [], [], []
    for _ in range(n):
        lat, lon = R.uniform(-89.9, 89.9), R.uniform(-180, 180)
        if R.random() < 0.2:
            lat = R.choice([R.uniform(86, 88), R.uniform(-0.2, 0.2), R.uniform(-88, -86)])
        surface = R.random() < 0.3
        odd = R.getrandbits(1)
        yz, xz = cpr_encode(lat, lon, odd, surface)
        if R.random() < 0.1:
            yz, xz = R.getrandbits(17), R.getrandbits(17)
        tc = R.choice([5, 6, 7, 8]) if surface else R.choice([9, 11, 18, 20])
        h = es(17, tc, me_pos(rbits(12), odd, yz, xz) if not surface else
               me_surface(R.randrange(128), 1, 0, odd, yz, xz))["hex"]
        span = 0.5 if surface else 2.0
        rlat = max(-90, min(90, lat + R.uniform(-span, span)))
        rlon = lon + R.uniform(-span, span)
        if R.random() < 0.1:
            rlat, rlon = float(R.randrange(-89, 90)), float(R.randrange(-180, 180))
        fn = pms.adsb.surface_position_with_ref if surface else pms.adsb.airborne_position_with_ref
        loc.append({"hex": h, "lat_ref": rlat, "lon_ref": rlon, "surface": surface, "pos": list(fn(h, rlat, rlon))})
    for _ in range(n):
        lat, lon = R.uniform(-89.9, 89.9), R.uniform(-180, 180)
        if R.random() < 0.2:
            lat = R.choice([R.uniform(86, 88), R.uniform(-0.2, 0.2), R.uniform(-88, -86), R.uniform(10.3, 10.6)])
        lat2, lon2 = lat + R.uniform(-0.05, 0.05), lon + R.uniform(-0.05, 0.05)
        e = cpr_encode(lat, lon, 0)
        o = cpr_encode(max(-89.99, min(89.99, lat2)), ((lon2 + 180) % 360) - 180, 1)
        if R.random() < 0.1:
            e = (R.getrandbits(17), R.getrandbits(17))
        he = es(17, 11, me_pos(rbits(12), 0, *e))["hex"]
        ho = es(17, 11, me_pos(rbits(12), 1, *o))["hex"]
        t0, t1 = R.choice([(1.0, 2.0), (2.0, 1.0), (1.5, 1.5)])
        a, b = (he, ho) if R.random() < 0.7 else (ho, he)
        ta, tb = (t0, t1) if a == he else (t1, t0)
        p = pms.adsb.airborne_position(a, b, ta, tb)
        glob.append({"a": a, "b": b, "ta": ta, "tb": tb, "pos": list(p) if p else None})
    edges = [0.0, 1e-9, 1e-7, -1e-9, 87.0, -87.0, 86.9995, 87.0005, 87.001, 86.535, 87.5, 90.0, -90.0, 45.0]
    for k in range(2, 59):                                       # NL transition latitudes
        a = 1 - math.cos(math.pi / 30)
        lt = math.degrees(math.acos(math.sqrt(a / (1 - math.cos(2 * math.pi / k)))))
        edges += [lt, math.nextafter(lt, 0), math.nextafter(lt, 90), -lt, lt + 1e-9, lt - 1e-9]
    edges += [R.uniform(-90, 90) for _ in range(300)]
    for x in edges:
        nl.append([x, int(pms.common.cprNL(x))])
    return {"local": loc, "global": glob, "nl": nl}


def main():
    ambiguous_pool()
    (OUT / "cpr.json").write_text(json.dumps(gen_cpr(800 * SCALE), separators=(",", ":")))
    fields = gen_fields(2400 * SCALE)
    (OUT / "fields.json").write_text(json.dumps(fields, separators=(",", ":")))
    cases = [scenario(k) for k in range(110 * SCALE)]
    (OUT / "tracks.json").write_text(json.dumps(cases, separators=(",", ":")))
    errs = sum(1 for c in cases if "error" in c)
    print(f"fields: {len(fields)} frames, {sum('field_error' in f for f in fields)} with field_error; "
          f"tracks: {len(cases)} cases, {errs} raised")
    (OUT / "demod.json").write_text(json.dumps(gen_demod(), separators=(",", ":")))
    # BDS50/60 resolution against an aircraft reference (pyModeS is50or60 through the python tracks)
    ambiguous_pool(60 + 200 * SCALE)
    res = []
    for mb in AMBIG:
        f = surv(R.choice([20, 21]), ac13(), rbits(14), mb)
        for _ in range(4):
            spd, trk_, alt = R.uniform(80, 520), R.uniform(0, 360), R.choice([R.uniform(0, 45000), 925, 37000])
            res.append({"hex": f["hex"], "spd": spd, "trk": trk_, "alt": alt,
                        "reg": pms.bds.is50or60(f["hex"], spd, trk_, alt)})
    (OUT / "is50or60.json").write_text(json.dumps(res, separators=(",", ":")))


if __name__ == "__main__":
    main()
