"""
dvbt2_p1.py -- DVB-T2 P1 preamble symbol (ETSI EN 302 755 clause 9.8; research/drone_transmitters/src_etsi/dvbt2.txt;
research, PLAN D5). Synthetic validation only (the spec's worked example and generated P1 symbols).

P1 = C | A | B: A is a 1K OFDM symbol (1024 T) with 384 active carriers (table 68) DBPSK-modulated by
{CSS_S1, CSS_S2, CSS_S1} (table 69) and scrambled (PRBS 1 + X^14 + X^15, init 100111001000110); C = the first 542 T of A
and B = the last 482 T of A, both shifted up by f_SH = 1 / (1024 T). T = 7 / (8 B) us (8 MHz: 7/64 us).
Detection: y[m] y*[m + 542] (C vs A) and y[m] y*[m - 482] (B vs A) rotate at f_SH; de-rotated and summed over the 542 / 482
samples they peak at the P1 start whatever the carrier offset. Decoding: FFT of A, integer carrier offset from the
active-carrier mask, descrambling, DBPSK soft bits between consecutive active carriers, correlation with the 8 S1 and 16 S2
patterns -> S1 (preamble format), S2 field 1 (FFT size / guard-interval class), S2 field 2 (mixed).
"""
from __future__ import annotations

import numpy as np

KP1 = np.array([
    44, 45, 47, 51, 54, 59, 62, 64, 65, 66, 70, 75, 78, 80, 81, 82, 84, 85, 87, 88, 89, 90, 94, 96, 97, 98, 102, 107,
    110, 112, 113, 114, 116, 117, 119, 120, 121, 122, 124, 125, 127, 131, 132, 133, 135, 136, 137, 138, 142, 144, 145,
    146, 148, 149, 151, 152, 153, 154, 158, 160, 161, 162, 166, 171,
    172, 173, 175, 179, 182, 187, 190, 192, 193, 194, 198, 203, 206, 208, 209, 210, 212, 213, 215, 216, 217, 218, 222,
    224, 225, 226, 230, 235, 238, 240, 241, 242, 244, 245, 247, 248, 249, 250, 252, 253, 255, 259, 260, 261, 263, 264,
    265, 266, 270, 272, 273, 274, 276, 277, 279, 280, 281, 282, 286, 288, 289, 290, 294, 299, 300, 301, 303, 307, 310,
    315, 318, 320, 321, 322, 326, 331, 334, 336, 337, 338, 340, 341, 343, 344, 345, 346, 350, 352, 353, 354, 358, 363,
    364, 365, 367, 371, 374, 379, 382, 384, 385, 386, 390, 395, 396, 397, 399, 403, 406, 411, 412, 413, 415, 419, 420,
    421, 423, 424, 425, 426, 428, 429, 431, 435, 438, 443, 446, 448, 449, 450, 454, 459, 462, 464, 465, 466, 468, 469,
    471, 472, 473, 474, 478, 480, 481, 482, 486, 491, 494, 496, 497, 498, 500, 501, 503, 504, 505, 506, 508, 509, 511,
    515, 516, 517, 519, 520, 521, 522, 526, 528, 529, 530, 532, 533, 535, 536, 537, 538, 542, 544, 545, 546, 550, 555,
    558, 560, 561, 562, 564, 565, 567, 568, 569, 570, 572, 573, 575, 579, 580, 581, 583, 584, 585, 586, 588, 589, 591,
    595, 598, 603, 604, 605, 607, 611, 612, 613, 615, 616, 617, 618, 622, 624, 625, 626, 628, 629, 631, 632, 633, 634,
    636, 637, 639, 643, 644, 645, 647, 648, 649, 650, 654, 656, 657, 658, 660, 661, 663, 664, 665, 666, 670, 672, 673,
    674, 678, 683,
    684, 689, 692, 696, 698, 699, 701, 702, 703, 704, 706, 707, 708, 712, 714, 715, 717, 718, 719, 720, 722, 723, 725,
    726, 727, 729, 733, 734, 735, 736, 738, 739, 740, 744, 746, 747, 748, 753, 756, 760, 762, 763, 765, 766, 767, 768,
    770, 771, 772, 776, 778, 779, 780, 785, 788, 792, 794, 795, 796, 801, 805, 806, 807, 809])
S1_HEX = ["124721741D482E7B", "47127421481D7B2E", "217412472E7B1D48", "742147127B2E481D",
          "1D482E7B12472174", "481D7B2E47127421", "2E7B1D4821741247", "7B2E481D74214712"]
S2_HEX = ["121D4748212E747B1D1248472E217B7412E247B721D174841DED48B82EDE7B8B",
          "4748121D747B212E48471D127B742E2147B712E2748421D148B81DED7B8B2EDE",
          "212E747B121D47482E217B741D12484721D1748412E247B72EDE7B8B1DED48B8",
          "747B212E4748121D7B742E2148471D12748421D147B712E27B8B2EDE48B81DED",
          "1D1248472E217B74121D4748212E747B1DED48B82EDE7B8B12E247B721D17484",
          "48471D127B742E214748121D747B212E48B81DED7B8B2EDE47B712E2748421D1",
          "2E217B741D124847212E747B121D47482EDE7B8B1DED48B821D1748412E247B7",
          "7B742E2148471D12747B212E4748121D7B8B2EDE48B81DED748421D147B712E2",
          "12E247B721D174841DED48B82EDE7B8B121D4748212E747B1D1248472E217B74",
          "47B712E2748421D148B81DED7B8B2EDE4748121D747B212E48471D127B742E21",
          "21D1748412E247B72EDE7B8B1DED48B8212E747B121D47482E217B741D124847",
          "748421D147B712E27B8B2EDE48B81DED747B212E4748121D7B742E2148471D12",
          "1DED48B82EDE7B8B12E247B721D174841D1248472E217B74121D4748212E747B",
          "48B81DED7B8B2EDE47B712E2748421D148471D127B742E214748121D747B212E",
          "2EDE7B8B1DED48B821D1748412E247B72E217B741D124847212E747B121D4748",
          "7B8B2EDE48B81DED748421D147B712E27B742E2148471D12747B212E4748121D"]
S1_NAMES = {0: "T2_SISO", 1: "T2_MISO", 2: "Non-T2 (FEF)", 3: "T2_LITE_SISO", 4: "T2_LITE_MISO"}
S2_BASE = {0: "2K, any GI", 1: "8K, GI 1/32 / 1/16 / 1/8 / 1/4", 2: "4K, any GI", 3: "1K, any GI", 4: "16K, any GI",
           5: "32K, GI 1/32 / 1/16 / 1/8", 6: "8K, GI 1/128 / 19/256 / 19/128", 7: "32K, GI 1/128 / 19/256 / 19/128"}
S2_LITE = {0: "2K, any GI", 1: "8K, GI 1/32 / 1/16 / 1/8 / 1/4", 2: "4K, any GI", 3: "16K, GI 1/128 / 19/256 / 19/128",
           4: "16K, GI 1/32 / 1/16 / 1/8 / 1/4", 6: "8K, GI 1/128 / 19/256 / 19/128"}
N1, LC, LB = 1024, 542, 482


def _hexbits(h):
    return np.array([(int(c, 16) >> (3 - i)) & 1 for c in h for i in range(4)], np.int8)


S1_BITS = np.array([_hexbits(h) for h in S1_HEX])          # 8 x 64
S2_BITS = np.array([_hexbits(h) for h in S2_HEX])          # 16 x 256


def scrambler(n=384, init="100111001000110"):
    """PRBS 1 + X^14 + X^15 (clause 5.2.4 generator) -> +-1 (0 -> +1)."""
    reg = [int(c) for c in init]
    out = np.empty(n, np.int8)
    for i in range(n):
        b = reg[13] ^ reg[14]
        out[i] = b
        reg = [b] + reg[:14]
    return 1 - 2 * out


SCR = scrambler()


def mss(s1, s2):
    """MSS_SCR (+-1, 384) for S1 (0..7) and S2 (0..15)."""
    seq = np.r_[S1_BITS[s1], S2_BITS[s2], S1_BITS[s1]]
    diff = np.cumprod(np.where(seq == 1, -1, 1))            # DBPSK from the dummy +1
    return diff * SCR


def p1_symbol(s1, s2):
    """Baseband P1 at the elementary rate 1/T (2048 samples, centred on DC)."""
    t = np.arange(N1)
    a = np.zeros(N1, complex)
    for k, c in zip(KP1, mss(s1, s2)):
        a += c * np.exp(2j * np.pi * (k - 426) * t / N1)
    a /= np.sqrt(384)
    sh = lambda x, t0: x * np.exp(2j * np.pi * (np.arange(len(x)) + t0) / N1)          # noqa: E731
    return np.r_[sh(a[:LC], 0), a, sh(a[N1 - LB:], LC + N1)]


def metric(y):
    """P1 correlation per start sample n (normalised, ~1 at a clean P1), over C-A and B-A."""
    y = np.asarray(y, complex)
    n = len(y) - 2 * N1
    if n <= 0:
        return np.zeros(0)
    m = np.arange(len(y))
    rot = np.exp(-2j * np.pi * m / N1)
    q1 = np.r_[y[:-LC] * np.conj(y[LC:]), np.zeros(LC)] * rot
    q2 = np.r_[np.zeros(LB), y[LB:] * np.conj(y[:-LB])] * rot
    e = np.abs(y) ** 2
    c1, c2, ce = (np.r_[0, np.cumsum(v)] for v in (q1, q2, e))
    s = np.arange(n)
    num = np.abs(c1[s + LC] - c1[s]) + np.abs(c2[s + 2 * N1] - c2[s + LC + N1])
    den = (ce[s + LC + N1] - ce[s]) * 0.5 + 1e-30                                          # ~ energy of 542 + 482 samples
    den = 0.5 * ((ce[s + LC] - ce[s]) + (ce[s + LC + N1] - ce[s + N1])
                 + (ce[s + 2 * N1] - ce[s + LC + N1]) + (ce[s + LC + N1 - LB] - ce[s + N1 - LB])) + 1e-30
    return num / den


def decode(y, s0):
    """P1 starting at sample s0 of y (elementary rate) -> dict(s1, s2, names, integer carrier offset, scores)."""
    a = np.asarray(y[s0 + LC:s0 + LC + N1], complex)
    if len(a) < N1:
        return None
    A = np.fft.fft(a)
    P = np.abs(A) ** 2
    best = None
    for sh in range(-100, 101):                               # +-100 carriers (+-0.9 MHz at 8 MHz)
        e = P[(KP1 - 426 + sh) % N1].sum()
        if best is None or e > best[0]:
            best = (e, sh)
    sh = best[1]
    r = A[(KP1 - 426 + sh) % N1] * SCR                        # descrambled: MSS_DIFF times the channel
    d = np.real(r[1:] * np.conj(r[:-1]))                      # + for bit 0, - for bit 1 (i = 1..383)
    soft = np.r_[0.0, d] / (np.abs(d).mean() + 1e-30)
    pm = lambda bits: 1.0 - 2.0 * bits                        # noqa: E731
    s1_sc = np.array([soft[1:64] @ pm(b)[1:] + soft[320:384] @ pm(b) for b in S1_BITS])
    s2_sc = np.array([soft[64:320] @ pm(b) for b in S2_BITS])
    s1, s2 = int(np.argmax(s1_sc)), int(np.argmax(s2_sc))
    srt1, srt2 = np.sort(s1_sc), np.sort(s2_sc)
    lite = s1 in (3, 4)
    out = dict(s1=f"{s1:03b}", s2=f"{s2:04b}", preamble=S1_NAMES.get(s1, "reserved"),
               fft_gi=(S2_LITE if lite else S2_BASE).get(s2 >> 1, "reserved") if s1 != 2 else
               ("undefined FEF part" if s2 >> 1 == 0 else "reserved"),
               mixed=bool(s2 & 1), carrier_offset=int(sh),
               s1_margin=round(float((srt1[-1] - srt1[-2]) / 127), 3), s2_margin=round(float((srt2[-1] - srt2[-2]) / 256), 3))
    return out


def find(y, thr=0.5):
    """P1 symbols in y (elementary rate): [(start, metric, decoded)]."""
    m = metric(y)
    out = []
    if not len(m):
        return out
    k = 0
    while k < len(m):
        j = int(np.argmax(m[k:k + 2 * N1])) + k
        if m[j] >= thr:
            d = decode(y, j)
            if d:
                out.append((j, float(m[j]), d))
            k = j + 2 * N1
        else:
            k += 2 * N1
    return out
