"""Bench (PLAN T2b, research/drone_transmitters/09_ble_coded_droneid_legacy.md): BLE LE Coded PHY S=8 synth +
decoder (coherent GFSK soft symbols -> S=8 pattern de-map -> K=4 Viterbi), synthetic only."""
import sys
import numpy as np
from spectrum_intel.receivers.rc_link import synth as S, gfsk_coherent as GC, detector as D

AA = 0x8E89BED6
PRE = [0, 0, 1, 1, 1, 1, 0, 0] * 10                    # 80 uncoded symbols
P0, P1 = [0, 0, 1, 1], [1, 1, 0, 0]                    # S=8 mapper


def conv(bits):
    s = [0, 0, 0]
    out = []
    for b in bits:
        out += [b ^ s[0] ^ s[1] ^ s[2], b ^ s[1] ^ s[2]]   # G0 = 1+x+x^2+x^3 (17), G1 = 1+x^2+x^3 (13)
        s = [b, s[0], s[1]]
    return out


def map8(cb):
    return [x for c in cb for x in (P1 if c else P0)]


def coded_bits(adv_addr, adv_data, ch=37, pdu_type=2):
    lsb = lambda by: [(v >> i) & 1 for v in by for i in range(8)]
    payload = bytes(adv_addr) + bytes(adv_data)
    hdr = bytes([pdu_type | 0x40, len(payload)])
    pdu = lsb(hdr) + lsb(payload)
    wh = S._ble_whiten(pdu + S._ble_crc_bits(pdu), ch)
    aa = [(AA >> i) & 1 for i in range(32)]
    b1 = map8(conv(aa + [0, 0] + [0, 0, 0]))           # AA + CI=00 (S=8) + TERM1
    b2 = map8(conv(list(wh) + [0, 0, 0]))              # PDU + CRC + TERM2
    return np.array(PRE + b1 + b2, np.uint8), len(pdu)


def viterbi(llr):
    """llr: (n, 2) soft values (+ -> 1) of coded bit pairs -> decoded bits (K=4, 8 states, terminated)."""
    n = len(llr)
    INF = 1e18
    pm = np.full(8, INF); pm[0] = 0
    prev = np.zeros((n, 8), np.int64); inb = np.zeros((n, 8), np.uint8)
    for t in range(n):
        npm = np.full(8, INF)
        for st in range(8):
            if pm[st] >= INF:
                continue
            s0, s1, s2 = (st >> 2) & 1, (st >> 1) & 1, st & 1
            for b in (0, 1):
                a0, a1 = b ^ s0 ^ s1 ^ s2, b ^ s1 ^ s2
                m = pm[st] - (2 * a0 - 1) * llr[t, 0] - (2 * a1 - 1) * llr[t, 1]
                ns = (b << 2) | (s0 << 1) | s1
                if m < npm[ns]:
                    npm[ns] = m; prev[t, ns] = st; inb[t, ns] = b
        pm = npm
    st = 0 if pm[0] < INF else int(np.argmin(pm))
    out = np.zeros(n, np.uint8)
    for t in range(n - 1, -1, -1):
        out[t] = inb[t, st]; st = prev[t, st]
    return out


def decode(y, s0, n_pdu_max=2 + 255):
    """y at 8 MSps, s0 = sample of preamble symbol 0."""
    nsym = 80 + 4 * 2 * 37 + 4 * 2 * (8 * n_pdu_max + 24 + 3)
    dm = GC.demodulate(y, s0, nsym)
    if dm is None:
        return None
    bits, rel = dm[0], dm[1]
    soft = (2.0 * bits - 1) * np.minimum(rel, np.r_[rel[1:], rel[-1]])     # symbol k: alpha_k from b_k, b_{k-1}
    soft = (2.0 * bits - 1) * np.minimum(rel, np.r_[rel[0], rel[:-1]])
    p = np.array(P1) * 2 - 1.0                          # +1 where the '1' pattern has a 1

    def llr_of(sym0, ncoded):
        s = soft[sym0:sym0 + 4 * ncoded].reshape(ncoded, 4) @ p
        return s.reshape(-1, 2)
    i = 80
    blk1 = viterbi(llr_of(i, 74))
    aa = sum(int(b) << k for k, b in enumerate(blk1[:32]))
    ci = int(blk1[32]) | (int(blk1[33]) << 1)
    i += 4 * 74
    avail = (len(soft) - i) // 8
    if avail < 16 + 24 + 3:
        return dict(aa=aa, ci=ci)
    blk2 = viterbi(llr_of(i, 2 * avail))
    d = D._ble_dewhiten(blk2, 37)
    ln = int(sum(int(d[8 + k]) << k for k in range(8)))
    n_pdu = 16 + 8 * ln
    ok = len(d) >= n_pdu + 24 and D._ble_crc_ok(d, n_pdu)
    return dict(aa=aa, ci=ci, len=ln, crc_ok=bool(ok), h=dm[2], pdu=D._ble_pdu(d, ln, 37, 37) if ok else None)


def burst(snr, rng, adv_addr, adv_data):
    sym, _ = coded_bits(adv_addr, adv_data)
    iq, _ = S.gfsk_burst(1e6, 250e3, bits=sym, bt=0.5, cfo_hz=rng.uniform(-40e3, 40e3), fs=8e6, rng=rng)
    x = np.r_[np.zeros(400), iq, np.zeros(400)]
    n0 = 1 / (1e6 * 10 ** (snr / 10))
    x = x + (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x))) * np.sqrt(n0 * 8e6 / 2)
    X = np.fft.fft(x); f = np.fft.fftfreq(len(x), 1 / 8e6); X[np.abs(f) > 0.75e6] = 0
    return np.fft.ifft(X), len(sym)


if __name__ == "__main__":
    sys.path.insert(0, "python/tests"); sys.path.insert(0, "python/tests/receivers")
    from test_gfsk_coherent import ADDR, location_advert
    adv = location_advert()
    for snr in [float(a) for a in sys.argv[1:]] or [10, 4, 0, -2, -4]:
        rng = np.random.default_rng(1)
        ok = 0
        for t in range(10):
            y, ns = burst(snr, rng, ADDR, adv)
            r = decode(y, 400, n_pdu_max=2 + len(adv) + 6)
            ok += bool(r and r.get("crc_ok"))
        print(f"SNR {snr:5.1f} dB  packets ok {ok}/10  (last: aa {r.get('aa', 0):08X} ci {r.get('ci')} h {r.get('h', 0):.2f})")
