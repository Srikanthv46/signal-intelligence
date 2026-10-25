"""Coherent BLE 1M GFSK detection + CRC-24 Chase correction (rc_link/gfsk_coherent.py; PLAN L1, bench/ble_lowsnr).
Synthetic only (research): constants regenerated from first principles, BER, realistic BLE 4 Remote ID adverts (one
25-byte message + counter per advert, as F3411 legacy advertising sends them) through rc_link detect at low SNR."""
import struct
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import erf

sys.path.insert(0, str(Path(__file__).resolve().parent / "receivers"))
from test_det_droneid import LAT, LON, UAS, _msg                          # noqa: E402

from spectrum_intel.receivers.rc_link import detector as D, gfsk_coherent as GC, synth as S   # noqa: E402

FC, FS = 2420e6, 61.44e6
AA = np.array([(D.BLE_ADV_AA >> i) & 1 for i in range(32)], np.uint8)
ADDR = bytes([0x11, 0x22, 0x33, 0x44, 0x55, 0xC6])


def c0_pulse(bt=0.5, sps=8, L=3):
    """Laurent's C0 for h 0.5 GFSK with a Gaussian frequency pulse truncated to L bits."""
    t = np.arange(0, (L + 1) * sps) / sps
    s = np.sqrt(np.log(2)) / (2 * np.pi * bt)
    tt = np.linspace(0, L, 4000)
    g = 0.5 * (erf((tt - L / 2 + 0.5) / (np.sqrt(2) * s)) - erf((tt - L / 2 - 0.5) / (np.sqrt(2) * s)))
    q = np.cumsum(g) * (tt[1] - tt[0]); q *= 0.5 / q[-1]
    qf = lambda x: np.interp(x, tt, q, left=0.0, right=0.5)          # noqa: E731

    def Sx(x):
        x = np.abs(x)
        out = np.zeros_like(x)
        m1 = x < L
        out[m1] = np.sin(np.pi * qf(x[m1])) / np.sin(np.pi / 2)
        m2 = (x >= L) & (x < 2 * L)
        out[m2] = np.sin(np.pi / 2 - np.pi * qf(x[m2] - L)) / np.sin(np.pi / 2)
        return out
    c = np.ones_like(t)
    for i in range(L):
        c *= Sx(t + i)
    return c


def test_constants():
    assert np.allclose(c0_pulse(), GC.C0, atol=1e-9)
    rng = np.random.default_rng(0)
    bits = rng.integers(0, 2, 4000).astype(np.uint8)
    iq, _ = S.gfsk_burst(1e6, 250e3, bits=bits, bt=0.5, fs=8e6, rng=rng)
    z = np.convolve(iq, GC.H_LPF, "same")
    d = np.angle(z[8:] * np.conj(z[:-8]))[np.arange(3999) * 8]
    a = 2 * bits.astype(int) - 1
    k = np.arange(1, 3997)
    for i in (0, 1):
        for j in (0, 1):
            kk = k[((a[k - 1] == a[k]) == bool(i)) & ((a[k + 1] == a[k]) == bool(j))]
            assert abs(np.mean(d[kk] * a[kk]) - GC.TAB[i, j]) < 1e-4          # empirical table (bench sample set)


def test_crc_syndrome():
    bits = np.array(S.ble_adv_bits(ADDR, b"\x02\x01\x06", ch=37), np.uint8)
    d = D._ble_dewhiten(bits[40:], 37)
    n_pdu = 16 + 8 * (6 + 3)
    assert D._ble_crc_ok(d, n_pdu)
    syn, contrib = GC.crc_syndrome_parts(d, n_pdu)
    assert syn == 0
    for i in (0, 17, n_pdu - 1, n_pdu, n_pdu + 23):
        e = d.copy(); e[i] ^= 1
        assert GC.crc_syndrome_parts(e, n_pdu)[0] == contrib[i] != 0
    rel = np.ones(40 + n_pdu + 24); rel[40 + 30] = 0.1; rel[40 + 50] = 0.2
    e = d.copy(); e[30] ^= 1; e[31] ^= 1; e[50] ^= 1; e[51] ^= 1      # two b-symbol errors
    r = GC.chase(e, n_pdu, rel, 40, protect=range(8, 16))
    assert r is not None and r[1] == 2 and np.array_equal(r[0], d)


def _burst(snr_db, rng, bits, dev=250e3, drift=0.0):
    iq, _ = S.gfsk_burst(1e6, dev, bits=bits, bt=0.5, cfo_hz=rng.uniform(-40e3, 40e3), fs=8e6, rng=rng)
    iq = iq * np.exp(1j * np.pi * drift * 1e6 * (np.arange(len(iq)) / 8e6) ** 2)
    x = np.r_[np.zeros(200), iq, np.zeros(200)]
    n0 = 1 / (1e6 * 10 ** (snr_db / 10))
    x = x + (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x))) * np.sqrt(n0 * 8e6 / 2)
    X = np.fft.fft(x); f = np.fft.fftfreq(len(x), 1 / 8e6); X[np.abs(f) > 0.75e6] = 0
    return np.fft.ifft(X)


def test_ber_8db():
    rng = np.random.default_rng(1)
    e = n = 0
    for _ in range(20):
        bits = np.r_[[1, 0] * 4, AA, rng.integers(0, 2, 360)].astype(np.uint8)
        y = _burst(8, rng, bits)
        i0 = GC.find_access_address(GC.diff_products(y), AA)[0][0]
        b = GC.demodulate(y, i0 - 64, len(bits))[0]
        m = min(len(b), len(bits)) - 2
        e += np.count_nonzero(b[40:m] != bits[40:m]); n += m - 40
    assert e / n < 2e-3, e / n


def basic_id_advert(counter=0):
    body = bytes([counter]) + _msg(0, 2, bytes([(1 << 4) | 2]) + UAS.encode() + b"\x00\x00\x00")
    return bytes([2, 0x01, 0x06]) + S.remote_id_adv_data(body)


def location_advert(counter=1):
    loc = _msg(1, 2, bytes([0x20, 90, 40]) + struct.pack("<b", 4)
               + struct.pack("<ii", int(LAT * 1e7), int(LON * 1e7))
               + struct.pack("<HHH", int(1120 / 0.5), int(1122 / 0.5), int(1118 / 0.5))
               + bytes([(5 << 4) | 11, (3 << 4) | 3]) + struct.pack("<H", 12345) + bytes([0, 0]))
    return bytes([2, 0x01, 0x06]) + S.remote_id_adv_data(bytes([counter]) + loc)


def capture(snr_db, seed=5, dev=250e3):
    rng = np.random.default_rng(seed)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k, adv in enumerate([basic_id_advert(0), location_advert(1), basic_id_advert(2)]):
        bits = np.array(S.ble_adv_bits(ADDR, adv, ch=37), np.uint8)
        iq, _ = S.gfsk_burst(1e6, dev, bits=bits, bt=0.5, cfo_hz=15e3, fs=FS, rng=rng)
        bursts.append(dict(t0_us=3000 + 9000 * k, f_hz=2402e6, iq=iq, bw_hz=1.0e6))
    return S.embed(x, FS, FC, bursts, snr_db)[0]


def _rid(em):
    return [e.decoded for e in em if (e.decoded or {}).get("ad", {}).get("remote_id")]


def test_realistic_adverts_low_snr():
    for snr, need in ((8, 3), (6, 2)):
        em = D.detect(capture(snr), FS, FC)
        got = _rid(em)
        assert len(got) >= need, (snr, len(got))
        for g in got:
            assert g["crc_ok"] and g["demod"] == "coherent" and abs(g["mod_index"] - 0.5) <= 0.03
            r = g["ad"]["remote_id"]
            assert r.get("uas_id", UAS) == UAS
            if "lat" in r:
                assert abs(r["lat"] - LAT) < 1e-5 and abs(r["lon"] - LON) < 1e-5


def test_off_nominal_index():
    em = D.detect(capture(8, seed=7, dev=228e3), FS, FC)          # h 0.456
    got = _rid(em)
    assert len(got) >= 2 and all(abs(g["mod_index"] - 0.456) <= 0.03 for g in got)   # h scatter ~0.02 on 376 bits


def test_no_false_decode():
    rng = np.random.default_rng(11)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k in range(3):                                            # 1 Mb/s GFSK without the advertising AA
        bits = rng.integers(0, 2, 400).astype(np.uint8)
        iq, _ = S.gfsk_burst(1e6, 250e3, bits=bits, bt=0.5, fs=FS, rng=rng)
        bursts.append(dict(t0_us=3000 + 9000 * k, f_hz=2402e6, iq=iq, bw_hz=1.0e6))
    y = S.embed(x, FS, FC, bursts, 8)[0]
    em = D.detect(y, FS, FC)
    assert not [e for e in em if (e.decoded or {}).get("crc_ok")]


def test_latency():
    rng = np.random.default_rng(3)
    bits = np.array(S.ble_adv_bits(ADDR, location_advert(), ch=37), np.uint8)
    y = _burst(8, rng, bits)
    D.ble_coherent_decode(y, 8e6, 37, AA)
    t = time.perf_counter()
    for _ in range(5):
        r = D.ble_coherent_decode(y, 8e6, 37, AA)
    dt = (time.perf_counter() - t) / 5
    assert r and r["crc_ok"] and dt < 0.02, dt
