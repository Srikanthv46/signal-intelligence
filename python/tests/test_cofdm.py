"""COFDM / DVB-T receiver (research: research/drone_transmitters/06_dvbt_cofdm.md) on synthetic DVB-T signals: numerology,
carrier, TPS decode (BCH) at several bandwidths / modes / guard intervals, down to 0 dB in-band SNR, and no false link on
noise or on a non-DVB-T OFDM signal of matching width."""
import numpy as np
import pytest

from spectrum_intel.receivers.cofdm import detector as D, synth as S, tables as T

FS, N = 61.44e6, 1 << 21


def make(B, mode, gi, qam=4, cr=0, snr=15.0, seed=1, off=5.3e6 + 1234, pilots=True):
    rng = np.random.default_rng(seed)
    L = T.MODES[mode]["N"] * (1 + 1 / T.GUARD[gi]) / (8e6 * B / 7)
    y = S.dvbt_capture_burst(int(N / FS / L) + 2, B, FS, mode=mode, gi=gi, qam=qam, start_symbol=17, rng=rng, cr_hp=cr,
                             constellation={4: 0, 16: 1, 64: 2}[qam], cell_id=0x1234, pilots=pilots)[:N]
    sig = np.zeros(N, np.complex64)
    sig[:len(y)] = y
    sig *= np.exp(2j * np.pi * off / FS * np.arange(N)).astype(np.complex64)
    nz = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * np.sqrt(FS / (B * 1e6) * 10 ** (-snr / 10) / 2)
    return (sig + nz).astype(np.complex64)


def test_tables_vectors():
    assert "".join(map(str, T.prbs(48))) == "111111111110000000001100000001111000001100110001"
    info = list(T.SYNC) + T.bits(0b010111, 6) + T.bits(0, 13) + T.bits(2, 2) + T.bits(0, 2) + T.bits(0, 14)
    assert "".join(map(str, T.bch_parity(info))) == "10101011100111"


@pytest.mark.parametrize("B,mode,gi,qam,cr,snr", [(8, "2K", "1/8", 4, 0, 15), (7, "2K", "1/32", 16, 4, 15),
                                                   (8, "2K", "1/8", 4, 0, 0)])
def test_tps_decoded(B, mode, gi, qam, cr, snr):
    em = D.detect(make(B, mode, gi, qam, cr, snr), FS, 900e6)
    assert len(em) == 1 and em[0].label == "dvbt"
    F = em[0].features
    assert F["mode"] == mode and F["guard_interval"] == gi and F["channel_bw_mhz"] == B
    assert abs(F["carrier_mhz"] - (900 + 5.3012340)) < 20e-6
    t = F["tps"]
    assert t["bch_ok"] and t["mode"] == mode and t["guard_interval"] == gi and t["code_rate_hp"] == T.CODE_RATE[cr]
    assert t["constellation"] == {4: "QPSK", 16: "16-QAM", 64: "64-QAM"}[qam]


@pytest.mark.parametrize("B,mode,gi", [(2, "2K", "1/16"), (4, "2K", "1/4"), (8, "8K", "1/32")])
def test_numerology_long_frames(B, mode, gi):
    em = D.detect(make(B, mode, gi), FS, 900e6)
    assert len(em) == 1 and em[0].label == "dvbt_like"            # the TPS frame is longer than the block
    F = em[0].features
    assert F["mode"] == mode and F["guard_interval"] == gi and F["channel_bw_mhz"] == B and "68 symbols" in F["tps_note"]


def test_no_false_link():
    rng = np.random.default_rng(3)
    noise = ((rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)).astype(np.complex64)
    assert D.detect(noise, FS, 900e6) == []
    # LTE-like OFDM, 15 kHz SCS, 300 subcarriers (4.5 MHz: the 5 MHz DVB-T width window), normal CP
    n_fft, cp = 4096, 288                                         # at 61.44 MSps
    sym = []
    for _ in range(N // (n_fft + cp)):
        X = np.zeros(n_fft, complex)
        X[np.r_[1:151, n_fft - 150:n_fft]] = (rng.choice([-1, 1], 300) + 1j * rng.choice([-1, 1], 300))
        s = np.fft.ifft(X)
        sym.append(np.r_[s[-cp:], s])
    y = np.concatenate(sym)
    y = y / np.sqrt(np.mean(np.abs(y) ** 2)) * 3
    x = noise.copy()
    x[:len(y)] += y.astype(np.complex64)
    assert D.detect(x, FS, 900e6) == []


def test_no_pilots_not_reported():
    # a CP-OFDM signal with the DVB-T 2K / 8 MHz numerology but no pilots (proprietary COFDM, or another periodic
    # signal folding at a DVB-T period: real-capture false positives on DJI OcuSync / 2.4 GHz ISM before the rule)
    assert D.detect(make(8, "2K", "1/8", pilots=False), FS, 900e6) == []


# ---- DVB-T2 P1 (research, PLAN D5; EN 302 755 clause 9.8)
from spectrum_intel.receivers.cofdm import dvbt2_p1 as P1   # noqa: E402


def test_p1_spec_example():
    m = P1.mss(0, 0)                                          # EN 302 755 9.8.2.2 example: S1 = 000, S2 = 0000
    assert list(P1.SCR[:4]) == [-1, 1, -1, 1] and list(m[:4]) == [-1, 1, -1, -1]
    assert m[63] == 1 and m[64] == -1 and m[320] == 1 and m[-1] == 1


@pytest.mark.parametrize("s1,s2", [(0, 0), (1, 3), (3, 9), (2, 0), (4, 12)])
@pytest.mark.parametrize("snr_db", [20, -3])
def test_p1_decode(s1, s2, snr_db):
    rng = np.random.default_rng(s1 * 16 + s2)
    p = P1.p1_symbol(s1, s2)
    y = np.r_[np.zeros(3000), p, np.zeros(3000)] * np.exp(2j * np.pi * 0.013 * np.arange(len(p) + 6000))
    y = y + (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y))) * np.sqrt(10 ** (-snr_db / 10) / 2 * 384 / 1024)
    f = P1.find(y, D.P1_THR)
    assert len(f) == 1 and abs(f[0][0] - 3000) <= 8                  # broad peak at -3 dB (0.9 us)
    d = f[0][2]
    assert d["s1"] == f"{s1:03b}" and d["s2"] == f"{s2:04b}" and d["carrier_offset"] == 13


def test_p1_no_false_alarm():
    rng = np.random.default_rng(2)
    assert P1.metric(rng.standard_normal(200000) + 1j * rng.standard_normal(200000)).max() < D.P1_THR / 2
    assert P1.metric(S.dvbt_baseband(40, mode="2K", gi="1/4", rng=rng)).max() < D.P1_THR / 2
