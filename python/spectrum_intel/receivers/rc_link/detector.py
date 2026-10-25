#!/usr/bin/env python3
"""
det_rc_link.py -- narrowband RC / telemetry control links (non-DJI drones).

    detect(x, fs, fc, ctx) -> list[Emission]       one Emission per burst
    group_hops(emissions)  -> list[dict]           RC-specific link grouping

Two independent detectors, one burst list:

  LoRa (ExpressLRS 2.4 GHz / 900 MHz, TBS Crossfire 50 Hz, any SX127x/SX128x CSS)
      Blind scan: the band is cut into sub-bands one LoRa BW wide (FFT-domain
      slice, critically sampled at fs = BW). At critical sampling a train of
      identical upchirps is exactly N-periodic, so multiplying by a periodic
      conjugate base chirp turns the whole preamble into ONE tone and a K-symbol
      FFT integrates it coherently (gain K*2^SF) regardless of symbol timing.
      Wi-Fi / OFDM / FSK / CW do not dechirp to a tone, so the test is robust in
      a busy band. Candidates are verified by finding the 2.25-downchirp SFD
      (opposite-polarity dechirp), which also separates timing from CFO:
          upchirp bin = cfo - tau,   downchirp bin = cfo + tau.
      Then the packet is walked symbol by symbol on the recovered grid:
      preamble length, sync-word symbols, payload length, SNR, and the drift of
      the fractional peak position -> chip-clock ppm.
  GFSK / FSK bursts (FrSky, FlySky, DSMX, ELRS FLRC/GFSK, Crossfire FSK, SiK)
      Matched-area energy detection on a normalised spectrogram at three
      bandwidth scales, narrowband test (flanks quiet while the burst is on),
      then per-burst channelisation: constant-envelope test (4th moment vs.
      the constant-envelope prediction at the measured SNR), instantaneous-
      frequency clusters (deviation, CFO), symbol rate from transition timing.
      1 Mb/s h~0.5 bursts are tried as BLE (preamble + access address
      0x8E89BED6, de-whitening, CRC-24); BLE/BR-like bursts are labelled
      bluetooth_like, never rc_fhss. 2.4 GHz bursts 0.5-3 MHz wide are sliced at
      2 Mb/s and need a >= 16-bit alternating preamble to be "ESB-like" links
      (esb.py); the packets of one link are CRC-16/CCITT-checked together.

Mirror images (IQ imbalance, 2*fc - f) of stronger bursts are dropped.

Fingerprints: cfo_ppm (carrier offset vs the ELRS channel grid / ppm; one
crystal -> one value across hops) is reliable; group_hops adds timer_ppm (packet
interval vs nominal = MCU timer crystal). clock_ppm is reported but biased (see
lora_emission).

Labels describe the physical layer; protocol names are hints ("_like") unless
decoded (BLE CRC) or matched on every ELRS parameter (SF, BW, preamble, length).
"""
from __future__ import annotations
import numpy as np
from scipy import ndimage, fft as sfft

from spectrum_intel.core import remote_id as RID
from . import bredr as BR
from . import sik as SIK
from . import gfsk_coherent as GC
from . import ble_coded as BC
from . import crossfire as CF
from . import elrs as ELRS
from . import fsk_frame as FF
from . import cc2500 as FRSKY
from . import afhds2a as AF2
from . import xn297 as XN
from . import radiolink as RL
from . import hubsan as HB
from . import hott as HOTT
from . import hitec as HITEC
from . import toys as TOYS
from . import ieee802154 as ZB
from . import flarm as FL
from . import adsl as ADSL
from . import ogn as OGN
from . import mlrs as ML
from spectrum_intel.core import bands_india as B
from spectrum_intel.core import common as C
from spectrum_intel.core.emission import Emission

from . import esb as ESB

# ============================================================== protocol tables
# ExpressLRS (src/src/common.cpp ExpressLRS_AirRateConfig, master Oct 2026):
# name: (band, sf, bw, cr_den, preamble, payload_bytes, hop_every, interval_us)
ELRS_LORA = {
    "2g4_500hz":     ("2G4", 5, 812.5e3, 6, 12, 8, 4, 2000),
    "2g4_333hz_8ch": ("2G4", 5, 812.5e3, 8, 12, 13, 4, 3003),
    "2g4_250hz":     ("2G4", 6, 812.5e3, 8, 14, 8, 4, 4000),
    "2g4_150hz":     ("2G4", 7, 812.5e3, 8, 12, 8, 4, 6666),
    "2g4_100hz_8ch": ("2G4", 7, 812.5e3, 8, 12, 13, 4, 10000),
    "2g4_50hz":      ("2G4", 8, 812.5e3, 8, 12, 8, 2, 20000),
    "900_200hz":     ("900", 6, 500e3, 7, 8, 8, 4, 5000),
    "900_100hz_8ch": ("900", 6, 500e3, 8, 8, 13, 4, 10000),
    "900_100hz":     ("900", 7, 500e3, 7, 8, 8, 4, 10000),
    "900_50hz":      ("900", 8, 500e3, 7, 10, 8, 4, 20000),
    "900_25hz":      ("900", 9, 500e3, 7, 10, 8, 2, 40000),
}
ELRS_DOMAINS = {   # FHSS.cpp: (f_start, f_stop, n_channels)
    "ISM2G4": (2400.4e6, 2479.4e6, 80), "FCC915": (903.5e6, 926.9e6, 40),
    "AU915": (915.5e6, 926.9e6, 20), "EU868": (863.275e6, 869.575e6, 13),
    "IN866": (865.375e6, 866.95e6, 4), "AU433": (433.42e6, 434.42e6, 3),
    "EU433": (433.1e6, 434.45e6, 3), "US433": (433.25e6, 438.0e6, 8),
    "US433W": (423.5e6, 438.0e6, 20), "TH920": (920.5e6, 924.7e6, 8),
}
# packet intervals of RC links we know (us) -- used by group_hops to name a rate
KNOWN_INTERVALS_US = {
    "elrs": [1000, 2000, 3003, 4000, 5000, 6666, 10000, 20000, 40000],
    "crossfire": [6666, 20000], "frsky_d16": [9000], "afhds2a": [3850], "dsmx": [11000, 22000],
}
# GFSK families: (label, symbol-rate range, deviation range) -- see report table
FSK_CLASSES = [
    ("elrs_flrc_like", (560e3, 740e3), (90e3, 260e3)),
    ("elrs_gfsk_like", (250e3, 350e3), (60e3, 140e3)),
    ("afhds2a_like", (420e3, 580e3), (80e3, 260e3)),
    ("crossfire_fsk_like", (78e3, 92e3), (30e3, 55e3)),
    ("frsky_d16_like", (62e3, 80e3), (40e3, 75e3)),
    ("sik_like", (40e3, 260e3), (15e3, 160e3)),
]
BLE_ADV_AA = 0x8E89BED6
# 2 Mb/s ESB-like packets: 2.4 GHz, a packet of 2-byte preamble + address + 1..255 bytes + CRC (~50 us .. 1.1 ms),
# 99 % bandwidth 0.5-3 MHz as measured (bench 2026-10-17: 0.53-2.1 MHz)
ESB_BAND = (2370e6, 2520e6)
ESB_US = (40.0, 1100.0)
ESB_BW = (0.5e6, 3.0e6)
ESB_MARGIN_US = 50.0
# GSM TDMA bursts (GMSK, constant envelope) look like a sub-GHz FSK packet: one time slot (or k adjacent slots) long,
# on the 200 kHz channel raster, inside a GSM-900 / 1800 band (traffic carriers have no FCCH for the GSM receiver)
GSM_SLOT_US = 576.9
GSM_SLOT_TOL_US = 60.0       # a single burst without its guard / ramps: 546-577 us
GSM_RASTER_TOL_HZ = 15e3     # receiver ppm (~10 ppm at 960 MHz) + carrier offset


def _lora_bws_for(fc):
    """Candidate (bw, SF list) per band: 2.4 GHz -> SX128x 812.5 kHz (ELRS) and
    1625 kHz; sub-GHz -> 500 kHz (ELRS 900, Crossfire LoRa) and 250/125 kHz, the
    125 kHz up to SF12 (LoRaWAN, Meshtastic, point-to-point sensors at 433/865 MHz).
    500 kHz SF5: SX126x telemetry links (live 890 MHz, testdata/live/LIVE_LORA_SF5_890MHz)."""
    if fc > 1.5e9:
        return [(812.5e3, [5, 6, 7, 8]), (1625e3, [6, 7, 8])]
    return [(500e3, [5, 6, 7, 8, 9]), (250e3, [7, 8, 9, 10]), (125e3, [7, 8, 9, 10, 11, 12])]


# ============================================================== small DSP helpers
def _slice(X, fs, f_off, half_bw, out_fs, t_adv=0.0):
    """FFT-domain channel slice centred on f_off. Signed bin indexing (the shared
    Channelizer.slice wraps negative offsets wrongly when L*out_fs/fs is not an
    integer). t_adv (output samples) applies a fractional time advance, y[n] =
    y(n + t_adv): critically-sampled LoRa payload symbols are only pure tones
    when sampled on the chip instants. Returns (y, exact output rate, exact f_off)."""
    L = len(X)
    M = int(round(L * out_fs / fs))
    k0 = int(round(f_off / fs * L))
    hb = int(half_bw / fs * L)
    k = np.arange(-hb, hb + 1)
    Y = np.zeros(M, np.complex64)
    v = X[(k0 + k) % L]
    if t_adv:
        v = v * np.exp(2j * np.pi * k / M * t_adv)
    Y[k % M] = v
    y = sfft.ifft(Y, workers=-1) * (M / L)
    return y.astype(np.complex64), M * fs / L, k0 * fs / L


def _upchirp(N):
    n = np.arange(N)
    return np.exp(1j * np.pi * (n * n / N - n)).astype(np.complex64)   # -BW/2 -> +BW/2, N-periodic


def _peak_interp(S, k):
    """Parabolic interpolation of |FFT|^2 peak (log domain) at integer bin k."""
    n = len(S)
    a, b, c = np.log(S[(k - 1) % n] + 1e-30), np.log(S[k] + 1e-30), np.log(S[(k + 1) % n] + 1e-30)
    d = a - 2 * b + c
    return k + (0.5 * (a - c) / d if d < 0 else 0.0)


def _wrap(v, period):
    return (v + period / 2) % period - period / 2


# ============================================================== LoRa: blind scan
# SF5 stage-1 threshold: N = 32 opens 4-16x more windows per second than SF7-9, so eta 24 passes ~10 noise
# peaks per 1.75 s and ~170 per 1.75 s of DMR/NXDN (all refine to nothing, 12 s). eta 40 is still ~-8 dB
# in-band SNR; SF5 decodes from ~0 dB (stage-1 ratio >= 108 at -2 dB, 614 on the live packet).
ETA_SF5 = 40.0


def lora_scan(X, fs, fc, bw, sfs, K=8, eta=24.0, chunk=48, edge=0.45, f_range=None):
    """Stage 1: coherent K-symbol dechirp FFT in every BW/4-spaced sub-band (only those
    centred inside f_range = (lo, hi) Hz offsets from fc, when given).
    Returns candidates dict(f_sb, sf, pol, win_start, ratio, bin)."""
    L = len(X)
    M = int(round(L * bw / fs)); out_fs = M * fs / L
    hb = int(bw / 2 / fs * L)
    centers = np.arange(-edge * fs + bw / 2, edge * fs - bw / 2 + 1, bw / 4)
    if f_range is not None:
        centers = centers[(centers >= f_range[0]) & (centers <= f_range[1])]
    cands = []
    refs = {sf: _upchirp(2 ** sf) for sf in sfs}
    for c0 in range(0, len(centers), chunk):
        cs = centers[c0:c0 + chunk]
        Y = np.zeros((len(cs), M), np.complex64)
        kk = np.arange(-hb, hb + 1)
        for i, f in enumerate(cs):
            k0 = int(round(f / fs * L))
            Y[i, kk % M] = X[(k0 + kk) % L]
        y = sfft.ifft(Y, axis=1, workers=-1)
        for sf in sfs:
            N = 2 ** sf
            nsym = M // N
            if nsym < K:                                  # capture shorter than K symbols
                continue
            u = np.tile(refs[sf], nsym)
            hop = N * (K // 2)                            # each sample in two windows
            for pol in (+1, -1):
                z = y[:, :nsym * N] * (np.conj(u) if pol > 0 else u)[None, :]
                W = np.lib.stride_tricks.sliding_window_view(z, K * N, axis=1)[:, ::hop, :]
                F = sfft.fft(W, axis=2, workers=-1)
                S = F.real ** 2 + F.imag ** 2
                # A fractional symbol-timing offset makes the de-chirped preamble an
                # N-periodic sequence whose energy sits on the two integer
                # (N-spaced) lines around -tau; CFO adds scalloping between FFT
                # bins. Tone energy = best 2x2 group {k, k+1, k+K, k+K+1}
                # (worst-case loss ~1.5 dB instead of ~7 dB); noise = rest.
                Sa = S + np.roll(S, -1, axis=2)
                S2 = Sa + np.roll(Sa, -K, axis=2)
                pk = S2.max(axis=2)
                r = pk * (K * N - 4) / (S.sum(axis=2) - pk + 1e-30)
                hit = np.argwhere(r > (max(eta, ETA_SF5) if sf == 5 else eta))
                for i, w in hit:
                    cands.append(dict(f_sb=float(cs[i]), sf=sf, bw=bw, pol=pol, win=int(w * hop),
                                      ratio=float(r[i, w]), bin=int(np.argmax(S2[i, w])), K=K))
    return cands, out_fs


def _nms_candidates(cands, bw):
    """Keep the strongest candidate per (time, freq) neighbourhood -- one LoRa
    packet lights up neighbouring sub-bands, SFs and windows."""
    cands = sorted(cands, key=lambda c: -c["ratio"])
    keep = []
    for c in cands:
        T = 2 ** c["sf"] / c["bw"]
        t0 = c["win"] / c["bw"]
        dup = False
        for k in keep:
            Tk = 2 ** k["sf"] / k["bw"]
            if abs(k["f_sb"] - c["f_sb"]) < 0.8 * max(bw, k["bw"]) and \
                    abs(k["win"] / k["bw"] - t0) < 16 * max(T, Tk):
                dup = True; break
        if not dup:
            keep.append(c)
    return keep


# ============================================================== LoRa: refine one packet
# ============================================================== LoRa: refine one packet
# Conventions: the reference chirp is N-periodic from stream sample 0, so after
# de-chirping, an upchirp whose symbol boundary is at sample b (tau = b mod N)
# is a tone at bin  cfo - pol*tau (+ pol*value for payload symbols) and a
# downchirp at cfo + pol*tau, independent of where the FFT window sits. Window
# alignment only matters for keeping one symbol per window.
def _sym_spec(z, start, N, pad=8):
    start = int(round(start))
    if start < 0 or start + N > len(z):
        return None
    seg = z[start:start + N]
    S = np.abs(sfft.fft(seg, N * pad)) ** 2
    E = float(np.sum(np.abs(seg) ** 2))
    k = int(np.argmax(S))
    pk = S[k]
    s2 = max((E - pk / N) / (N - 1), 1e-30)          # noise / sample (signal energy removed)
    return S, _peak_interp(S, k) / pad, pk / (N * s2), s2   # spectrum, fine peak (bins), peak/noise-bin, s2


def _streams(X, fs, f_off, bw, N, pol, t_adv=0.0, cache=None):
    """cache: dict t_adv -> _slice output for this f_off/bw (candidates of one sub-band share it: the
    full-capture slice is most of a refinement's cost, and it does not depend on SF or polarity)."""
    if cache is None:
        y, ofs, f_act = _slice(X, fs, f_off, bw / 2, bw, t_adv)
    else:
        if t_adv not in cache:
            cache[t_adv] = _slice(X, fs, f_off, bw / 2, bw, t_adv)
        y, ofs, f_act = cache[t_adv]
    uu = np.resize(_upchirp(N), len(y))
    return (y * (np.conj(uu) if pol > 0 else uu), y * (uu if pol > 0 else np.conj(uu)), ofs, f_act)


def lora_refine(X, fs, fc, cand, n_iter=3, min_sfd_ratio=5.0, cache=None):
    """Stage 2: SFD search, CFO/timing split, aligned refinement, symbol walk."""
    sf, bw, pol, K = cand["sf"], cand["bw"], cand["pol"], cand["K"]
    N = 2 ** sf
    z_up, z_dn, ofs, f_act = _streams(X, fs, cand["f_sb"], bw, N, pol, cache=cache)
    n = len(z_up)
    s0 = cand["win"]
    S = np.abs(sfft.fft(z_up[s0:s0 + K * N], K * N * 8)) ** 2
    k_u = (_peak_interp(S, int(np.argmax(S))) / (8 * K)) % N
    # preamble extent: the de-chirped preamble tone (bin k_u) is present in any
    # N-window inside the preamble, whatever its alignment
    def r_at(z, st, kexp, tol=2):
        p = _sym_spec(z, st, N)
        if p is None:
            return 0.0
        kb = int(round(kexp * 8)) % (8 * N)
        return float(p[0][(kb + np.arange(-tol, tol + 1)) % (8 * N)].max() / (N * p[3]))
    prof = {j: r_at(z_up, s0 + j * N, k_u) for j in range(-20, 34)}
    # a real preamble puts the same tone in every symbol of the stage-1 window;
    # interference that only looks coherent over K symbols (a transient, a
    # burst edge) does not -> cheap early reject before the expensive steps
    # the stage-1 window may straddle the preamble start or end (NMS keeps one window per
    # packet): test the K-symbol run with the most preamble in it, up to K/2 symbols away
    j0 = max(range(-(K // 2), K // 2 + 1), key=lambda j: (sum(prof[j + i] >= 3.0 for i in range(K)), -abs(j)))
    rin = np.array([prof[j0 + i] for i in range(K)])
    if np.sum(rin >= 3.0) < K - 2 or rin.max() > 6 * max(np.median(rin), 1e-9):
        return None
    j_end = j0 + K - 1
    while j_end + 2 in prof and (prof[j_end + 1] + prof[j_end + 2]) / 2 >= 3.0:
        j_end += 1
    pre_end = s0 + (j_end + 1) * N                       # +-1 symbol
    best = (0.0, None, None)
    for adv in (0.0, 0.5):                               # half-sample grid: SFD lines split less
        zd = z_dn if adv == 0 else _streams(X, fs, cand["f_sb"], bw, N, pol, t_adv=adv, cache=cache)[1]
        for p in range(int(pre_end + 0.5 * N), int(pre_end + 3.75 * N), max(N // 4, 1)):
            if p + 2 * N > n:
                break
            Sd = np.abs(sfft.fft(zd[p:p + 2 * N], 2 * N * 4)) ** 2
            r = Sd.max() / (np.median(Sd) / np.log(2))
            if r > best[0]:
                best = (r, p + adv, Sd)
    r_d, p_d, Sd = best
    j_start = j0
    while j_start - 2 in prof and (prof[j_start - 1] + prof[j_start - 2]) / 2 >= 3.0:
        j_start -= 1
    partial = dict(partial=True, sf=sf, bw=bw, pol=pol, N=N, ofs=ofs, f_sb=fc + cand["f_sb"], k_u=k_u,
                   pre_start=s0 + j_start * N, pre_end=pre_end, n_pre_est=j_end - j_start + 1,
                   ratio1=cand["ratio"], K=K, sfd_window_in_capture=int(pre_end + 3.75 * N) + 2 * N <= n)
    if p_d is None or r_d < min_sfd_ratio:
        return partial
    k_d = (_peak_interp(Sd, int(np.argmax(Sd))) / 8) % N
    cfo0 = _wrap((k_u + k_d) / 2, N / 2)
    # cfo is ambiguous by N/2 (tau with it); the SFD grid point by +-1 symbol.
    # Score each hypothesis: downchirp power at the predicted bin in the two SFD
    # symbols, minus the symbols either side (must NOT be full downchirps),
    # plus preamble power 3..5 symbols earlier.
    hyps = []
    for cfo in (cfo0, cfo0 + N / 2 if cfo0 < 0 else cfo0 - N / 2):
        tau = (pol * (k_d - cfo)) % N
        zu2, zd2, _, f2 = _streams(X, fs, f_act + cfo * ofs / N, bw, N, pol)
        c2 = cfo - (f2 - f_act) / (ofs / N)
        kdx, kux = c2 + pol * tau, c2 - pol * tau
        j0 = np.floor((p_d - tau) / N)
        best = None
        for j in (j0 - 1, j0, j0 + 1):
            b = tau + j * N
            sc = (r_at(zd2, b, kdx, 1) + r_at(zd2, b + N, kdx, 1) - r_at(zd2, b - N, kdx, 1)
                  - r_at(zd2, b + 2 * N, kdx, 1) + sum(r_at(zu2, b - m * N, kux, 1) for m in (3, 4, 5)))
            if best is None or sc > best[0]:
                best = (sc, cfo, b)
        hyps.append(best)
    # the score above usually settles the N/2 ambiguity, but not when stage 1 kept the
    # sub-band BW/2 off the carrier (a ghost of the same chirp): then the wrong pair also
    # walks as a "packet". Walk both and keep the one with the stronger preamble + SFD.
    walks, failed_walk = [], False
    for sc, cfo, b_sfd in sorted(hyps, key=lambda h: -h[0]):
        wk = _refine_walk(X, fs, fc, f_act, ofs, cfo, b_sfd, N, bw, pol, n_iter)
        if wk:
            walks.append(wk)
        elif wk is False:
            failed_walk = True
    if not walks:
        return partial if failed_walk else None
    walk = max(walks, key=lambda wk: wk["z"])
    walk.update(sfd_ratio=r_d, ratio1=cand["ratio"], prem=partial)
    return walk


def _refine_walk(X, fs, fc, f_act, ofs, cfo, b_sfd, N, bw, pol, n_iter):
    """Iterative CFO / SFD-timing refinement for one hypothesis, then the symbol walk.
    None: the refinement lost the packet; False: refined but the walk found no packet."""
    sf = int(np.log2(N))
    f_ref = f_act
    for it in range(n_iter):
        f_new = f_ref + cfo * ofs / N                     # re-centre the slice on the carrier
        z_up, z_dn, ofs, f_act = _streams(X, fs, f_new, bw, N, pol)
        cfo -= (f_act - f_ref) / (ofs / N); f_ref = f_act
        pu = [_sym_spec(z_up, b_sfd - m * N, N) for m in range(3, 9)]
        pdn = [_sym_spec(z_dn, b_sfd + m * N, N) for m in range(2)]
        pu = [p for p in pu if p is not None]; pdn = [p for p in pdn if p is not None]
        if len(pu) < 3 or not pdn:
            return None
        ku = _circ_mean([p[1] for p in pu], N, [p[2] for p in pu])
        kd = _circ_mean([p[1] for p in pdn], N, [p[2] for p in pdn])
        tb = b_sfd % N
        du = _wrap(ku - (cfo - pol * tb), N)
        dd = _wrap(kd - (cfo + pol * tb), N)
        cfo += (du + dd) / 2
        b_sfd += pol * (dd - du) / 2
    # final pass on the chip-aligned stream (fractional timing removed)
    frac = b_sfd - np.floor(b_sfd)
    z_up, z_dn, ofs, f_act = _streams(X, fs, f_ref, bw, N, pol, t_adv=frac)
    walk = _walk_packet(z_up, z_dn, np.floor(b_sfd), N, cfo, pol)
    if walk is None:
        return False
    walk["start_samp"] += frac; walk["end_samp"] += frac
    walk.update(f_carrier=fc + f_act + cfo * ofs / N, ofs=ofs, sf=sf, bw=bw, pol=pol, N=N, partial=False)
    return walk


def _circ_mean(v, period, w=None):
    v = np.asarray(v, float); w = np.ones_like(v) if w is None else np.asarray(w, float)
    a = np.angle(np.sum(w * np.exp(2j * np.pi * v / period)))
    return (a / (2 * np.pi) * period) % period


def _walk_packet(z_up, z_dn, b_sfd, N, cfo, pol, max_pre=40, max_pay=400):
    """Walk the symbol grid fixed by the SFD start b_sfd (fractional sample)."""
    base = cfo - pol * (b_sfd % N)                       # preamble tone bin
    syms = []                                            # (g, fractional deviation, peak/noise)
    pre, miss = [], 0
    for m in range(3, max_pre + 3):                      # preamble, walking back past 2 sync symbols
        p = _sym_spec(z_up, b_sfd - m * N, N)
        if p is None:
            break
        S, pk, r, s2 = p
        kb = int(round(base * 8)) % (8 * N)
        r_exp = S[(kb + np.arange(-1, 2)) % (8 * N)].max() / (N * s2)   # power at the expected bin
        if r_exp < 4.0:
            miss += 1
            if miss > 1:
                break
            continue
        miss = 0
        dv = _wrap(pk - base, N)
        pre.append((m, dv if abs(dv) < 0.5 else 0.0, r_exp, abs(dv) < 0.5 and r > 6))
    if len(pre) < 2:
        return None
    n_pre = max(q[0] for q in pre) - 2
    syms += [(-m, dv, r) for m, dv, r, good in pre if good]
    sfd_r = []
    for m in range(2):
        p = _sym_spec(z_dn, b_sfd + m * N, N)
        if p is not None:
            kb = int(round((cfo + pol * (b_sfd % N)) * 8)) % (8 * N)
            sfd_r.append(p[0][(kb + np.arange(-1, 2)) % (8 * N)].max() / (N * p[3]))
    sync_vals = []
    for m in (2, 1):
        p = _sym_spec(z_up, b_sfd - m * N, N)
        if p is not None:
            sync_vals.append(int(round(pol * _wrap(p[1] - base, N))) % N)
    p0 = b_sfd + 2.25 * N
    # the payload grid starts 2.25 symbols after the SFD, a quarter symbol off the
    # preamble grid, so its zero-value tone sits N/4 away from the preamble tone
    base_pay = base - pol * N / 4
    pay = []
    for j in range(max_pay):
        p = _sym_spec(z_up, p0 + j * N, N)
        if p is None:
            break
        pay.append((j, p[1], p[2]))
    if not pay:
        return None
    r = np.array([q[2] for q in pay])
    # payload end: per-symbol peak ratio is ~ln(N)+2 on noise (measured) and
    # ~1+N*snr on signal (snr from the preamble, 15% purity loss); stop where
    # the 4-symbol mean falls near noise, then CUSUM at the midpoint
    c_noise = np.log(N) + 2.0
    r_pre0 = np.array([q[2] for q in pre])
    sig = max(0.85 * np.mean(r_pre0), c_noise + 1.0)
    c = c_noise + 0.5 * (sig - c_noise)
    ma = np.convolve(np.r_[r, np.full(3, c_noise)], np.ones(4) / 4, "valid")
    run_end = next((i for i in range(len(ma)) if ma[i] < c_noise + 0.25 * (sig - c_noise)), len(r))
    run_end = min(run_end + 3, len(r))
    # a following packet on the same channel (ELRS: 2-symbol gap at 150 Hz)
    # starts with identical preamble symbols; random payload almost never
    # repeats one value (+-1) `run` times: (3/N)^(run-1) per symbol, so small N
    # needs a longer run (SF5: a 3-run cut half of all 85-symbol packets short)
    vals = np.array([int(round(pol * _wrap(q[1] - base_pay, N))) % N for q in pay])
    run = 3 if N >= 128 else 4 if N >= 64 else 5
    for j in range(len(vals) - run + 1):
        if all(abs(_wrap(vals[j + i] - vals[j], N)) <= 1 for i in range(1, run)) and min(r[j:j + run]) > c:
            run_end = min(run_end, j)
            break
    cs = np.cumsum(r[:max(run_end, 1)] - c)
    n_pay = int(np.argmax(cs)) + 1
    pay_vals = []
    for j, pk, rr in pay[:n_pay]:
        v = pol * _wrap(pk - base_pay, N)
        vi = np.round(v)
        pay_vals.append(int(vi) % N)
        syms.append((2.25 + j, pol * (v - vi), rr))
    # in-band SNR from the preamble power at the known bin (unbiased by argmax):
    # r = 1 + N*snr
    r_pre = np.array([q[2] for q in pre])
    snr_db = 10 * np.log10(max((np.mean(r_pre) - 1) / N, 1e-4))
    # verification z-score: preamble + SFD power at the predicted bins vs exp(1) noise
    n_ev = len(r_pre) + len(sfd_r)
    z_pre = float((r_pre.sum() + sum(sfd_r) - n_ev) / np.sqrt(n_ev))
    z_sfd = float((sum(sfd_r) - len(sfd_r)) / np.sqrt(max(len(sfd_r), 1)))
    # chip-clock drift: fractional peak deviation vs symbol index (weighted LS)
    g = np.array([q[0] for q in syms], float); d = np.array([q[1] for q in syms]); w = np.array([q[2] for q in syms])
    ok = w > 8
    slope = slope_err = np.nan
    if ok.sum() >= 6 and np.ptp(g[ok]) > 8:
        G = np.c_[np.ones(ok.sum()), g[ok]]
        W = w[ok] / w[ok].max()
        A = G.T @ (G * W[:, None])
        coef = np.linalg.solve(A, G.T @ (W * d[ok]))
        res = d[ok] - G @ coef
        s2 = np.sum(W * res ** 2) / max(ok.sum() - 2, 1)
        slope, slope_err = coef[1], float(np.sqrt(np.linalg.inv(A)[1, 1] * s2))
    return dict(start_samp=b_sfd - (n_pre + 2) * N, end_samp=p0 + n_pay * N, n_pre=n_pre, n_pay=n_pay,
                sync_vals=sync_vals, pay_vals=pay_vals, snr_db=snr_db, slope=slope, slope_err=slope_err,
                pre_ratio=float(np.median(r_pre)), pay_ratio=float(np.median(r[:n_pay])), z=z_pre, z_sfd=z_sfd,
                sfd_r=sfd_r, all_vals=vals.tolist())


def lora_payload_symbols(sf, pl_bytes, cr_den, crc=False, implicit=True, ldro=False):
    """Semtech SX127x time-on-air payload symbol count (SX128x SF5/6 differ slightly; unverified)."""
    num = 8 * pl_bytes - 4 * sf + 28 + 16 * crc - 20 * implicit
    return 8 + max(int(np.ceil(num / (4 * (sf - 2 * ldro)))) * cr_den, 0)


def _mlrs_match(sf, bw, n_pre, n_pay, fc_hz):
    """mLRS air modes consistent with the measured LoRa parameters (preamble 12, 91-byte frames; research, PLAN R3).
    The payload count uses the SX127x formula; the SX126x / SX1280 framings add a few symbols, hence the wider
    tolerance."""
    out = []
    band = "2G4" if fc_hz > 1.5e9 else "900"
    for name, (b, msf, mbw, cr, pre, pl, ival, dec) in ML.MLRS_LORA.items():
        if b != band or msf != sf or abs(mbw - bw) > 1e3:
            continue
        exp_pay = lora_payload_symbols(sf, pl, cr)
        if abs(n_pre - pre) <= 3 and abs(n_pay - exp_pay) <= 0.1 * exp_pay + 4:
            out.append(dict(mode=name, interval_us=ival, decodable=dec, exp_payload_syms=exp_pay))
    return out


def _elrs_match(sf, bw, n_pre, n_pay, fc_hz):
    """ELRS air-rate modes consistent with the measured LoRa parameters."""
    out = []
    for name, (band, msf, mbw, cr, pre, pl, hop, ival) in ELRS_LORA.items():
        if msf != sf or abs(mbw - bw) > 1e3:
            continue
        # SX127x counts programmed preamble + 4.25 internally -- accept +-3 symbols
        exp_pay = lora_payload_symbols(sf, pl, cr)
        out.append(dict(mode=name, interval_us=ival, hop_every=hop, preamble=pre, exp_payload_syms=exp_pay,
                        pre_ok=abs(n_pre - pre) <= 3, pay_ok=abs(n_pay - exp_pay) <= 4))
    return out


def _elrs_ota(w, sf, bw, matches):
    from . import lora_phy as P
    for cr_den, pl in dict.fromkeys((ELRS_LORA[m["mode"]][3], ELRS_LORA[m["mode"]][5]) for m in matches):
        d = P.decode_implicit(w["all_vals"], sf, cr_den - 4, pl, P.default_ldro(sf, bw))
        if not d:
            continue
        e = ELRS.parse(d["payload"])
        if e and e.get("crc_ok"):
            e.update(cr=d["cr"], ldro=False, sync_word=None, implicit_header=True, fec_corrected=d["fec_corrected"],
                     fec_detected_errors=d["fec_detected_errors"])
            return e
    return None


def _nominal_channel(f_hz):
    """Nearest ELRS channel (any domain) -> (domain, index, f_nominal)."""
    best = None
    for dom, (f0, f1, n) in ELRS_DOMAINS.items():
        sp = (f1 - f0) / (n - 1)
        i = int(round((f_hz - f0) / sp))
        if 0 <= i < n:
            fn = f0 + i * sp
            if best is None or abs(f_hz - fn) < abs(f_hz - best[2]):
                best = (dom, i, fn)
    return best


def lora_emission(w, fc, fs):
    sf, bw, N, ofs = w["sf"], w["bw"], w["N"], w["ofs"]
    f = w["f_carrier"]
    # chip clock: TX symbol = N(1+delta) of our samples; up-chirp peaks drift by -N*delta per symbol.
    # CAVEAT (synthetic validation): this estimate is biased (~+10 ppm) and
    # compressed (gain ~0.3-0.4) and its scatter is 10-50 ppm -- not usable as a
    # fingerprint yet; use cfo_ppm (per burst) and group_hops' timer_ppm (per link).
    clock_ppm, clock_err = None, None
    if np.isfinite(w["slope"]):
        delta = -w["slope"] / N * (1 if w["pol"] > 0 else -1)
        eps = ofs / (bw * (1 + delta)) - 1
        clock_ppm = eps * 1e6
        clock_err = w["slope_err"] / N * 1e6
    dec = _lora_payload(w, sf, bw)
    if dec and (dec.get("header_ok") or dec.get("implicit_header")):
        # the explicit header (or a CRC-verified implicit length) fixes the packet length:
        # trust it over the energy-based end
        w = dict(w, n_pay=dec["n_symbols"], end_samp=w["end_samp"] - w["n_pay"] * N + dec["n_symbols"] * N)
    nom = _nominal_channel(f)
    matches = [] if dec else _elrs_match(sf, bw, w["n_pre"], w["n_pay"], f)   # ELRS uses implicit headers
    full = [m for m in matches if m["pre_ok"] and m["pay_ok"]]
    if not dec and matches and f < 1.5e9:
        # ExpressLRS sub-GHz (SX127x, standard coding): implicit-header payload of each matching mode; a SYNC packet
        # verifies itself (its own UID4 / UID5 seed the CRC) -> hop index, nonce, rate, TLM ratio, UID bytes
        # (research, PLAN R2c). 2.4 GHz uses the SX1280's undocumented long-interleaving code: not decoded.
        dec = _elrs_ota(w, sf, bw, matches)
    # mLRS (research, PLAN R3 / T6): never IQ-inverted, preamble 12, 91-byte frames; the SX127x-framed SF6 mode decodes
    ml = [] if (dec or full or w["pol"] < 0) else _mlrs_match(sf, bw, w["n_pre"], w["n_pay"], f)
    if ml and any(m["decodable"] for m in ml):
        from . import lora_phy as P
        d = P.decode_implicit(w["all_vals"], sf, 1, ML.FRAME_LEN, False)
        dec = ML.parse(d["payload"]) if d else None
        if dec:
            dec.update(cr=d["cr"], ldro=False, sync_word=None, implicit_header=True, fec_corrected=d["fec_corrected"])
    if ml:
        family, label = "rc_fhss", f"mlrs{'' if (dec or {}).get('crc_ok') else '_like'}_lora_sf{sf}_bw{bw / 1e3:g}k"
        cfo = None; conf = 0.95 if (dec or {}).get("crc_ok") else 0.7
    elif nom is not None and abs(f - nom[2]) < 150e3 and full:
        family, label = "rc_fhss", f"elrs_lora_sf{sf}_bw{bw / 1e3:g}k"
        cfo = f - nom[2]; conf = 0.97 if (dec or {}).get("crc_ok") else 0.9
    elif nom is not None and abs(f - nom[2]) < 150e3 and matches:
        family, label = "rc_fhss", f"elrs_like_lora_sf{sf}_bw{bw / 1e3:g}k"
        cfo = f - nom[2]; conf = 0.7
    else:
        family, label = "lora", f"lora_sf{sf}_bw{bw / 1e3:g}k"
        cfo = None; conf = 0.6
        if dec and dec.get("crc_ok"):
            conf = 0.97
            if dec.get("protocol"):
                label = f"{dec['protocol'].lower()}_{label}"
        elif dec:
            conf = 0.8                                     # header checksum passed, payload unverified
    t0 = w["start_samp"] / ofs * 1e6; t1 = w["end_samp"] / ofs * 1e6
    feats = dict(modulation="lora", sf=sf, bw_hz=bw, bw_mhz=bw / 1e6, iq_inverted=w["pol"] < 0,
                 preamble_syms=w["n_pre"], payload_syms=w["n_pay"], sync_syms=w["sync_vals"],
                 symbol_us=N / bw * 1e6, burst_us=t1 - t0, f_carrier_hz=f,
                 cfo_hz=cfo, cfo_ppm=(cfo / f * 1e6) if cfo is not None else None,
                 clock_ppm=clock_ppm, clock_ppm_err=clock_err,
                 elrs_channel=(nom[0], nom[1]) if (nom and cfo is not None) else None,
                 elrs_modes=[m["mode"] for m in (full or matches)], mlrs_modes=[m["mode"] for m in ml],
                 preamble_ratio=w["pre_ratio"], sfd_ratio=w["sfd_ratio"], z_pre=w["z"], z_sfd=w["z_sfd"])
    if dec:
        feats.update(sync_word=dec.get("sync_word"), coding_rate=dec["cr"], ldro=dec["ldro"])
    return Emission(family=family, label=label, start_us=t0, end_us=t1,
                    f_lo_mhz=(f - bw / 2) / 1e6, f_hi_mhz=(f + bw / 2) / 1e6, center_mhz=f / 1e6,
                    snr_db=float(w["snr_db"]), confidence=conf, features=feats, decoded=dec)


def _lora_payload(w, sf, bw):
    """Explicit-header payload decode of a walked packet (lora_phy); SF5: implicit header,
    CRC-verified. None when the header block does not check (implicit-header link such as
    ExpressLRS, or noise)."""
    from . import lora_phy as P
    d = P.decode_sf5_implicit(w["all_vals"], w["n_pay"]) if sf == 5 else P.decode(w["all_vals"], sf, bw)
    if d is None:
        return None
    # a 5-bit header checksum passes on ~1 in 40 implicit-header packets by chance:
    # without a CRC pass, the header's length must also match the measured burst
    if not d.get("crc_ok") and abs(d["n_symbols"] - w["n_pay"]) > 4:
        return None
    sw = P.sync_word(w["sync_vals"])
    d["sync_word"] = f"0x{sw:02X}" if sw is not None else None
    d["network"] = P.SYNC_WORDS.get(sw)
    if d.get("complete"):
        pay = bytes.fromhex(d["payload_hex"])
        d["payload_ascii"] = P.ascii_view(pay)
        if d.get("crc_ok"):
            hdr = P.identify(sw, pay)
            if hdr:
                d["protocol"] = hdr.pop("protocol")
                d["fields"] = hdr
        if "protocol" not in d:
            # MAVLink in the payload of a transparent LoRa serial link (Ebyte E22 / E32 / E220 and similar;
            # research, PLAN G2): self-verifying through the MAVLink CRC_EXTRA, with or without the LoRa CRC
            mv = P.mavlink_in(pay)
            if mv:
                d["protocol"] = "MAVLink over LoRa (transparent serial link)"
                d["fields"] = mv
    return d


def lora_preamble_emission(p, fc):
    """Preamble-only detection (SFD not confirmed, e.g. below the SFD's own
    sensitivity): carrier known only to the sub-band (+-BW/4), no CFO split."""
    sf, bw, N, ofs = p["sf"], p["bw"], p["N"], p["ofs"]
    t0 = p["pre_start"] / ofs * 1e6; t1 = p["pre_end"] / ofs * 1e6
    f = p["f_sb"]
    nom = _nominal_channel(f)
    elrs = nom is not None and abs(f - nom[2]) < bw / 4 and any(
        m[1] == sf and abs(m[2] - bw) < 1e3 for m in ELRS_LORA.values())
    feats = dict(modulation="lora", sf=sf, bw_hz=bw, bw_mhz=bw / 1e6, iq_inverted=p["pol"] < 0,
                 preamble_syms=p["n_pre_est"], symbol_us=N / bw * 1e6, burst_us=t1 - t0, sfd_found=False,
                 cfo_hz=None, cfo_ppm=None, clock_ppm=None, f_uncertainty_hz=bw / 4,
                 stage1_ratio=p["ratio1"], elrs_channel=(nom[0], nom[1]) if elrs else None)
    return Emission(family="rc_fhss" if elrs else "lora",
                    label=f"{'elrs_like_' if elrs else ''}lora_preamble_sf{sf}_bw{bw / 1e3:g}k",
                    start_us=t0, end_us=t1, f_lo_mhz=(f - bw / 2) / 1e6, f_hi_mhz=(f + bw / 2) / 1e6,
                    center_mhz=f / 1e6, snr_db=None, confidence=0.5, features=feats)


def detect_lora(X, fs, fc, bws=None, eta=24.0, z_min=12.0, zsfd_min=4.0, eta_pre=60.0, npre_range=(5, 40),
                log=None, debug=None, f_range=None):
    """Blind LoRa detection over the band. Full detections need preamble + SFD +
    payload evidence; preamble-only detections need a much stronger stage-1
    coherent peak (eta_pre, K=8 symbols), a plausible preamble run length, and an SFD
    that is either past the end of the capture or seen (only the payload cut short).
    Sync symbols, when read, must be multiples of 8 (+-1): a sync word's nibbles x 8.
    debug (list) collects every refined candidate's statistics for calibration."""
    out = []
    for bw, sfs in (bws or _lora_bws_for(fc)):
        cands, _ = lora_scan(X, fs, fc, bw, sfs, eta=eta, f_range=f_range)
        cands = _nms_candidates(cands, bw)
        if log: log(f"  lora bw {bw/1e3:g}k: {len(cands)} stage-1 candidates")
        # refine sub-band by sub-band, sharing its slices (results kept in candidate order)
        walks, caches = [None] * len(cands), {}
        for i in sorted(range(len(cands)), key=lambda i: cands[i]["f_sb"]):
            f = cands[i]["f_sb"]
            caches = caches if f in caches else {f: {}}         # one sub-band's slices in memory at a time
            walks[i] = lora_refine(X, fs, fc, cands[i], cache=caches[f])
        for c, w in zip(cands, walks):
            if w is None:
                continue
            p = w if w["partial"] else w["prem"]
            # a sync word's two symbols are its nibbles x 8 (SX127x / SX126x / gr-lora_sdr): anything else is not LoRa
            sync_ok = w["partial"] or all(min(v % 8, 8 - v % 8) <= 1 for v in (w.get("sync_vals") or []))
            ok = (not w["partial"]) and w["z"] >= z_min and w["z_sfd"] >= zsfd_min and \
                w["n_pre"] >= 4 and w["n_pay"] >= 4 and sync_ok
            # preamble-only: a real preamble is followed by its SFD, so either the SFD lies past the end of the
            # capture, or it was seen (z_sfd) and only the payload is cut short
            sfd_ok = not p.get("sfd_window_in_capture", True) or (not w["partial"] and w["z_sfd"] >= zsfd_min)
            ok_pre = (not ok) and c["ratio"] >= eta_pre and npre_range[0] <= p["n_pre_est"] <= npre_range[1] \
                and sfd_ok and sync_ok
            if debug is not None:
                debug.append(dict(ok=ok, ok_pre=ok_pre, ratio1=c["ratio"], n_pre_est=p["n_pre_est"],
                                  sfd_window_in_capture=p.get("sfd_window_in_capture"),
                                  z=w.get("z"), z_sfd=w.get("z_sfd"), n_pre=w.get("n_pre"), n_pay=w.get("n_pay"),
                                  sf=c["sf"], bw=bw, f=w.get("f_carrier", p["f_sb"]),
                                  t=w.get("start_samp", p["pre_start"]) / p["ofs"] * 1e6))
            if ok:
                out.append(lora_emission(w, fc, fs))
            elif ok_pre:
                out.append(lora_preamble_emission(p, fc))
    return _dedup(out)


def _dedup(ems, f_tol_frac=0.5):
    # a decoded packet beats a copy of itself walked one symbol off (sync read as [0, 8])
    ems = sorted(ems, key=lambda e: (-bool((e.decoded or {}).get("crc_ok")), -bool(e.decoded),
                                     -(e.features.get("preamble_ratio") or 0)))
    keep = []
    for e in ems:
        bw = (e.f_hi_mhz - e.f_lo_mhz)
        # a chirp shifted in frequency is the same chirp shifted in time, so a confirmed
        # packet's preamble also lights up the neighbouring sub-bands: a preamble-only
        # detection within one BW of it is that ghost
        tol = 1.0 if e.features.get("sfd_found") is False else f_tol_frac
        if any(abs(k.center_mhz - e.center_mhz) < tol * bw and
               min(k.end_us, e.end_us) - max(k.start_us, e.start_us) > 0.3 * (e.end_us - e.start_us)
               for k in keep):
            continue
        keep.append(e)
    return sorted(keep, key=lambda e: e.start_us)


# ============================================================== FSK / GFSK bursts
def fsk_scan(x, fs, fc, nfft=512, k_sigma=8.0, widths=(1, 3, 7), t_frames=12, max_bw=3.0e6,
             min_us=30.0, edge=0.45):
    """Narrowband burst candidates from a normalised spectrogram (power over the
    per-bin floor), box-averaged at several bandwidths: a burst of bandwidth B
    and duration T is found at the SNR where B*T cells of excess integrate
    above k_sigma, i.e. matched-area energy detection."""
    hop = nfft // 2
    P, t_us, f_mhz = C.spectrogram(x, fs, fc, nfft, hop)
    fl = np.percentile(P, 10, axis=0) / 0.10536          # 10th pct of exp(1) = 0.105 x mean
    fl = ndimage.median_filter(fl, 15)
    Z = P / fl[None, :]
    keep = np.abs(f_mhz - fc / 1e6) <= edge * fs / 1e6
    # second normalisation: power over the local (+-3 MHz, per frame) median, so
    # a narrowband burst riding on a wideband signal (Wi-Fi, OcuSync video) still
    # stands out; the flank test below keeps wideband edges out
    loc = ndimage.median_filter(P, (1, 49), mode="nearest")
    Zl = P / np.maximum(loc, fl[None, :])
    E = np.full(Z.shape, -np.inf, np.float32)
    for ZZ, tfs in ((Z, (t_frames, 4 * t_frames, 16 * t_frames)), (Zl, (t_frames, 4 * t_frames))):
        for wf in widths:
            for tf in tfs:                                # 50, 200 (, 800) us integration
                Zf = ndimage.uniform_filter(ZZ, (tf, wf), mode="nearest")
                neff = wf * tf / 1.5                      # Hann 50% overlap -> correlated cells
                E = np.maximum(E, ((Zf - 1) * np.sqrt(neff)).astype(np.float32))
    E[:, ~keep] = -np.inf
    dt = hop / fs * 1e6; df = fs / nfft / 1e6
    cands = []
    claimed = np.zeros(E.shape, bool)
    # ascending thresholds: a very strong burst is merged with its own spectral
    # splatter / ADC-compression stripes at k_sigma; at 10x / 100x the core
    # separates out
    gap = max(1, int(round(60.0 / (hop / fs * 1e6))))   # bridge <=60 us dips (long low-SNR bursts fragment)
    for thr in (k_sigma, 10 * k_sigma, 100 * k_sigma, 1000 * k_sigma):
      mk = ndimage.binary_closing(E > thr, np.ones((gap, 1))) & np.isfinite(E)
      lab, n = ndimage.label(mk, np.ones((3, 3)))
      for i, sl in enumerate(ndimage.find_objects(lab)):
        ti, fi = sl
        if claimed[sl][lab[sl] == i + 1].any():
            continue
        dur = (ti.stop - ti.start) * dt
        bw = (fi.stop - fi.start) * df
        if dur < min_us or bw > max_bw / 1e6 or dur > 0.9 * t_us[-1]:
            continue
        sub = Z[sl] * (lab[sl] == i + 1)
        spec = sub.sum(0)
        fc_b = float(np.sum(f_mhz[fi] * spec) / max(spec.sum(), 1e-12))
        # flanks while the burst is on: a narrowband burst leaves them near the floor
        # flanks while the burst is on: in-band power must stand >= 2x (3 dB-ish
        # in excess terms) above BOTH flanks -- true for a narrowband burst on
        # the floor or on top of a wideband signal, false for pieces of one
        zin = Z[ti, fi].mean()
        lo = (f_mhz >= f_mhz[fi.start] - 2.0) & (f_mhz < f_mhz[fi.start] - 0.4) & keep
        hi = (f_mhz > f_mhz[fi.stop - 1] + 0.4) & (f_mhz <= f_mhz[fi.stop - 1] + 2.0) & keep
        fl_z = [np.median(Z[ti][:, m].mean(0)) for m in (lo, hi) if m.sum() >= 3]
        if not fl_z or (max(fl_z) - 1 > 0.3 * (zin - 1) and zin < 2.0 * max(fl_z)):
            continue
        claimed[sl] |= lab[sl] == i + 1
        cands.append(dict(t0_us=float(t_us[ti.start] - dt), t1_us=float(t_us[ti.stop - 1] + dt),
                          f_lo=float(f_mhz[fi.start] - df / 2) * 1e6, f_hi=float(f_mhz[fi.stop - 1] + df / 2) * 1e6,
                          fc_b=fc_b * 1e6, zin=float(zin), peak_e=float(E[sl].max())))
    return cands


def _kmeans2(v, it=20):
    a, b = np.percentile(v, 25), np.percentile(v, 75)
    for _ in range(it):
        m = np.abs(v - a) < np.abs(v - b)
        if m.all() or (~m).all():
            break
        a, b = v[m].mean(), v[~m].mean()
    return min(a, b), max(a, b)


def _symbol_rate(ifr, ofs, r_lo=20e3, r_hi=1.3e6):
    """Symbol rate of a binary FSK discriminator signal: spectral line of the
    transition sequence, refined by the timing-phase concentration of the
    transition instants. Returns (rate, concentration, phase)."""
    best = None
    for sm in (max(1, int(ofs / 4e6)), max(1, int(ofs / 1e6)), max(1, int(ofs / 250e3))):
        v = ndimage.uniform_filter1d(ifr, sm)
        bsg = np.sign(v - np.median(v))
        tr = np.flatnonzero(np.diff(bsg) != 0)
        if len(tr) < 8:
            continue
        ind = np.zeros(len(v)); ind[tr] = 1
        F = np.abs(np.fft.rfft(ind - ind.mean(), 4 * len(ind))) ** 2
        ff = np.fft.rfftfreq(4 * len(ind), 1 / ofs)
        sel = (ff > r_lo) & (ff < min(r_hi, ofs / (2.5 * sm)))
        if not sel.any():
            continue
        k = np.argmax(F[sel]); R = ff[sel][k]
        t = tr / ofs
        grid = R * (1 + np.linspace(-0.03, 0.03, 241))
        conc = np.abs(np.exp(2j * np.pi * np.outer(grid, t)).mean(1))
        j = int(np.argmax(conc))
        if best is None or conc[j] > best[1]:
            ph = np.angle(np.exp(2j * np.pi * grid[j] * t).mean())
            best = (float(grid[j]), float(conc[j]), float(ph))
    return best if best else (None, 0.0, 0.0)


# ---- BLE (Core spec Vol 6 Part B): LSB-first bits, whitening x^7+x^4+1, CRC-24
def _ble_dewhiten(bits, ch):
    reg = [1] + [(ch >> (5 - i)) & 1 for i in range(6)]
    out = np.empty(len(bits), np.uint8)
    for i, b in enumerate(bits):
        o = reg[6]
        out[i] = b ^ o
        reg = [o] + reg[:6]
        reg[4] ^= o
    return out


def _ble_crc_ok(bits, n_pdu_bits, init=0x555555):
    """Bit-serial CRC-24 (poly 0x00065B, LFSR per spec, init 0x555555) over the PDU;
    compares with the 24 transmitted CRC bits."""
    reg = [(init >> i) & 1 for i in range(24)]            # position 0 = LSB of init (spec loading)
    for b in bits[:n_pdu_bits]:
        fb = reg[23] ^ b
        reg = [fb] + reg[:23]
        for p in (1, 3, 4, 6, 9, 10):
            reg[p] ^= fb
    crc_tx = list(bits[n_pdu_bits:n_pdu_bits + 24])
    return len(crc_tx) == 24 and crc_tx == reg[::-1]        # transmitted position 23 first


def _ble_channel_index(f_hz):
    k = (f_hz - 2402e6) / 2e6
    ki = int(round(k))
    if abs(k - ki) > 0.1 or not 0 <= ki <= 39:
        return None
    if ki == 0: return 37
    if ki == 12: return 38
    if ki == 39: return 39
    return ki - 1 if ki < 12 else ki - 2


def _ble_pdu(d, ln, ch, wch):
    """Dewhitened PDU bits d (CRC verified) -> the decoded dict."""
    ptype = int(sum(int(d[i]) << i for i in range(4)))
    by = [int(sum(int(d[16 + 8 * j + i]) << i for i in range(8))) for j in range(ln)]
    out = dict(aa=f"{BLE_ADV_AA:08X}", crc_ok=True, channel_index_rf=ch, whitening_channel=wch,
               channel_mismatch=(wch != ch), pdu_type=ptype, pdu_len=ln,
               pdu_name={0: "ADV_IND", 1: "ADV_DIRECT_IND", 2: "ADV_NONCONN_IND", 3: "SCAN_REQ",
                         4: "SCAN_RSP", 5: "CONNECT_IND", 6: "ADV_SCAN_IND",
                         7: "ADV_EXT_IND/AUX_*"}.get(ptype, "?"),
               tx_addr_random=bool(d[6]))
    if ptype in (0, 2, 4, 6) and ln >= 6:
        out["adv_addr"] = ":".join(f"{v:02X}" for v in by[5::-1])
        out["adv_data_hex"] = bytes(by[6:]).hex()
        ad = RID.parse_ble_ad(bytes(by[6:]))       # AD structures; Remote ID (F3411) if present
        if ad:
            out["ad"] = ad
    elif ptype == 7 and ln >= 1:
        # extended advertising (BT5: ADV_EXT_IND on 37-39 pointing to an AUX_ADV_IND that carries the AdvData;
        # Remote ID over BT5 Long Range); research, PLAN T2b
        ext = BC.parse_extended(by)
        out["pdu_name"] = "ADV_EXT_IND" if ch in (37, 38, 39) else "AUX_ADV_IND"
        data = ext.pop("adv_data", b"")
        if ext.get("adv_addr"):
            out["adv_addr"] = ext.pop("adv_addr")
        out["extended"] = ext
        if data:
            out["adv_data_hex"] = data.hex()
            ad = RID.parse_ble_ad(data)
            if ad:
                if ad.get("remote_id"):
                    ad["remote_id"]["transport"] = "Bluetooth LE extended advertising (BT5)"
                out["ad"] = ad
    return out


def ble_try_decode(y, ofs, f_hz, t_rel=None):
    """Demodulate 1 Mb/s GFSK, look for preamble + advertising access address
    0x8E89BED6, de-whiten with the channel index of f_hz, check CRC-24.
    At 8 samples/bit the coherent detector with CRC Chase correction runs first (gfsk_coherent: ~4-5 dB more
    sensitive, ~5-10 ms; whitening channel of f_hz only); the hard-decision slicer is the fallback (any of the 40
    whitening channels: shifted RF labels, ~100 ms in Python).
    Returns dict or None."""
    ch = _ble_channel_index(f_hz)
    aa_bits = np.array([(BLE_ADV_AA >> i) & 1 for i in range(32)], np.uint8)
    coh = ble_coherent_decode(y, ofs, ch, aa_bits) or ble_coded_decode(y, ofs, ch, aa_bits)
    if coh:
        return coh
    sps = ofs / 1e6
    ifr = np.angle(y[1:] * np.conj(y[:-1]))
    pw = np.abs(y[1:]) ** 2
    on = pw > 0.25 * np.percentile(pw, 90)
    fallback = None
    for sm, ph in [(sm, ph) for sm in (max(1, int(sps / 2)), max(1, int(sps * 0.7))) for ph in np.arange(0, sps, max(sps / 8, 1))]:
        ifs = ndimage.uniform_filter1d(ifr, sm)
        thr = np.mean(_kmeans2(ifs[on])) if on.sum() > 16 else 0.0
        idx = (ph + np.arange(int((len(ifs) - ph) / sps)) * sps).astype(int)
        bits = (ifs[idx] > thr).astype(np.uint8)
        if len(bits) < 80:
            return fallback
        s = bits.astype(np.int8) * 2 - 1
        a = aa_bits.astype(np.int8) * 2 - 1
        corr = np.correlate(s, a, "valid")
        pos = np.flatnonzero(corr >= 30)                    # <=1 bit error in 32
        for p in pos[:3]:
            pdu = bits[p + 32:]
            if len(pdu) < 16 + 24:
                continue
            # whitening seed = channel index of f_hz; if the CRC fails try all
            # 40 (a pass with another index means the RF axis or the emitter's
            # channel map is shifted -- reported, not hidden)
            order = ([ch] if ch is not None else []) + [c for c in range(40) if c != ch]
            for wch in order:
                d = _ble_dewhiten(pdu, wch)
                ln = int(sum(int(d[8 + i]) << i for i in range(8)))
                n_pdu = 16 + 8 * ln
                if len(d) < n_pdu + 24 or not _ble_crc_ok(d, n_pdu):
                    continue
                return _ble_pdu(d, ln, ch, wch)
        if len(pos) and fallback is None:
            fallback = dict(aa=f"{BLE_ADV_AA:08X}", crc_ok=False, channel_index_rf=ch)
    return fallback


def ble_coded_decode(y, ofs, ch, aa_bits):
    """LE Coded PHY (Long Range, BT5 Remote ID; research, PLAN T2b): coded preamble + access address by soft
    correlation, coherent soft symbols, S=8 / S=2 de-mapping, soft Viterbi, de-whitening, CRC-24 (no Chase: Viterbi
    errors come in bursts). Returns the decoded dict (phy "LE Coded S=8" / "S=2") or None."""
    if ch is None or abs(ofs / (GC.SPS * 1e6) - 1) > 1e-3:
        return None
    p = GC.diff_products(y)
    for i0, nc in BC.find_coded(p, aa_bits):
        r = BC.decode(y, i0, aa_bits)
        if not r or r["aa_bit_errors"] > 4 or "bits" not in r:
            continue
        d = _ble_dewhiten(r["bits"], ch)
        if len(d) < 16 + 24:
            continue
        ln = int(sum(int(d[8 + i]) << i for i in range(8)))
        n_pdu = 16 + 8 * ln
        if len(d) < n_pdu + 24 or not _ble_crc_ok(d, n_pdu):
            continue
        out = _ble_pdu(d, ln, ch, ch)
        out.update(phy="LE Coded S=8" if r["ci"] == 0 else "LE Coded S=2", demod="coherent",
                   aa_bit_errors=r["aa_bit_errors"], mod_index=round(r["h"], 2), aa_ncorr=round(nc, 3))
        return out
    return None


def ble_coherent_decode(y, ofs, ch, aa_bits):
    """Low-SNR path (PLAN L1): access address by soft correlation, coherent detection (gfsk_coherent), CRC with
    Chase correction (<= 3 b-symbol flips among the 16 least reliable, length field protected) on the whitening
    channel of the RF frequency. Returns the decoded dict (demod "coherent", bit_corrections, mod_index) or None."""
    if ch is None or abs(ofs / (GC.SPS * 1e6) - 1) > 1e-3:          # the slice rate is 8 MSps to FFT-bin rounding
        return None
    p = GC.diff_products(y)
    for i0, nc in GC.find_access_address(p, aa_bits):
        nb_pre = min(8, i0 // GC.SPS)
        s0 = i0 - nb_pre * GC.SPS
        res = GC.demodulate(y, s0, nb_pre + 32 + 16 + 8 * 255 + 24)
        if res is None:
            continue
        bits, rel, h, cfo, dt = res
        # the access address on the coherent bit grid (+-1 bit: the soft search and the timing search may differ)
        sh = min((-1, 0, 1), key=lambda q: np.count_nonzero(bits[nb_pre + q:nb_pre + q + 32] != aa_bits)
                 if nb_pre + q >= 0 else 99)
        aa_err = int(np.count_nonzero(bits[nb_pre + sh:nb_pre + sh + 32] != aa_bits))
        if aa_err > 4:
            continue
        off = nb_pre + sh + 32
        pdu = bits[off:]
        if len(pdu) < 16 + 24:
            continue
        d = _ble_dewhiten(pdu, ch)
        ln = int(sum(int(d[8 + i]) << i for i in range(8)))
        n_pdu = 16 + 8 * ln
        if len(d) < n_pdu + 24:
            continue
        r = GC.chase(d, n_pdu, rel, off, protect=range(8, 16))
        if r is None:
            continue
        d2, nfix = r
        out = _ble_pdu(d2, ln, ch, ch)
        out.update(demod="coherent", bit_corrections=int(nfix), aa_bit_errors=aa_err, mod_index=round(h, 2),
                   aa_ncorr=round(nc, 3))
        return out
    return None


def fsk_refine(X, fs, fc, cand, noise_win_us=2000.0):
    """Channelise one candidate and measure it. Returns a list (a candidate
    box can hold several same-channel packets separated by short gaps)."""
    bw0 = cand["f_hi"] - cand["f_lo"]
    half = max(0.75 * bw0 + 60e3, 150e3)
    # never above the capture rate: on a 2.5 MSps capture the 2-16 MHz slice was an
    # up-sampling (12.8 M-point IFFTs per candidate, minutes per snapshot)
    ofs = float(min(max(8 * half, 2e6), 16e6, fs))
    y, ofs, f_act = _slice(X, fs, cand["fc_b"] - fc, half, ofs)
    a, b = int(cand["t0_us"] * 1e-6 * ofs), int(cand["t1_us"] * 1e-6 * ofs)
    a, b = max(a, 0), min(b, len(y))
    if b - a < 16:
        return []
    # pass 1: occupied bandwidth of the burst (PSD minus the local noise PSD),
    # then re-slice tightly so a 70 kBd burst is not measured in a 1 MHz slice
    nseg = 256
    def psd(v):
        m = len(v) // nseg
        if m < 1:
            return None
        return np.mean(np.abs(np.fft.fftshift(np.fft.fft(v[:m * nseg].reshape(m, nseg) * np.hanning(nseg),
                                                          axis=1), axes=1)) ** 2, 0)
    w0 = int(noise_win_us * 1e-6 * ofs)
    ps, pn = psd(y[a:b]), psd(np.r_[y[max(0, a - w0):max(0, a - 32)], y[min(len(y), b + 32):min(len(y), b + w0)]])
    if ps is not None and pn is not None:
        ex = ndimage.uniform_filter1d(np.clip(ps - pn, 0, None), 5)
        frq = np.fft.fftshift(np.fft.fftfreq(nseg, 1 / ofs))
        sel = np.flatnonzero(ex > 0.1 * ex.max())
        if len(sel):
            f_lo_b, f_hi_b = frq[sel[0]], frq[sel[-1]]
            half = max(0.6 * (f_hi_b - f_lo_b) + 2 * ofs / nseg, 60e3)
            ofs2 = float(min(max(8 * half, 1e6), 16e6, fs))
            y, ofs2, f_act = _slice(X, fs, f_act + (f_lo_b + f_hi_b) / 2, half, ofs2)
            a, b = int(a / ofs * ofs2), int(b / ofs * ofs2)
            ofs = ofs2
    pw = np.abs(y) ** 2
    w = int(noise_win_us * 1e-6 * ofs)
    out_idx = np.r_[max(0, a - w):max(0, a - int(10e-6 * ofs)), min(len(y), b + int(10e-6 * ofs)):min(len(y), b + w)]
    n0 = float(np.median(pw[out_idx]) / np.log(2)) if len(out_idx) > 32 else float(np.percentile(pw, 20) / 0.223)
    # precise edges: smoothed power crossing halfway (in dB-linear terms) between
    # noise and the burst level
    k = max(1, int(10e-6 * ofs))
    env = np.convolve(pw, np.ones(k) / k, "same")
    lvl = np.percentile(env[a:b], 75)
    if lvl < 1.15 * n0:
        # a short strong packet in a long box: fsk_scan's 800 us integration smears a 1 kHz train of ~70-100 us
        # packets (2 Mb/s ESB-like links) into ~870 us boxes, where it fills < 25 % of the box. Its level is then
        # the top of the box; a smoothed noise-only box stays below 1.3 n0 there, so 4 n0 (6 dB) keeps it out.
        lvl = np.percentile(env[a:b], 98)
        if lvl < 4 * n0:
            return []
    thr = n0 + 0.5 * (lvl - n0)
    base = max(0, a - 5 * k)
    on = env[base:min(len(y), b + 5 * k)] > thr
    on = ndimage.binary_closing(on, np.ones(max(1, int(8e-6 * ofs))))   # bridge <8 us dips
    lab_r, nr = ndimage.label(on)
    runs = [(base + r.start, base + r.stop) for (r,) in ndimage.find_objects(lab_r)
            if (r.stop - r.start) / ofs * 1e6 >= 20]
    out = [m for m in (_fsk_measure(X, fs, fc, y, ofs, f_act, half, n0, i0, i1) for i0, i1 in runs) if m]
    f_box = fc + f_act
    ble_box = 2370e6 <= f_box <= 2520e6 and not any((m.get("ble") or {}).get("crc_ok") for m in out)
    cf_box = 850e6 <= f_box <= 935e6 and not any(m.get("crossfire") for m in out)
    if len(runs) >= 2 and (b - a) / ofs >= 400e-6 and (ble_box or cf_box):
        # a long packet splits into pieces: at the half-level edges near the noise (BT5 Long Range at ~1 dB: ~70
        # pieces of an 8 ms AUX_ADV_IND) or where an S=2 block follows the S=8 one (500 kBd -> 1 MBd apparent rate):
        # one decode attempt over the whole box (research, PLAN T2b)
        m = _fsk_measure(X, fs, fc, y, ofs, f_act, half, n0, max(a, 0), min(b, len(y)))
        if m and ((m.get("ble") or {}).get("crc_ok") or m.get("crossfire")):
            return [m]
    return out


def _fsk_measure(X, fs, fc, y, ofs, f_act, half, n0, i0, i1):
    seg = y[i0:i1]
    if len(seg) < 16:
        return None
    m2 = float(np.mean(np.abs(seg) ** 2)); s_pow = max(m2 - n0, 1e-30)
    snr_slice = s_pow / n0
    # constant-envelope test: E|y|^4 for A*e^{j phi} + CN(0,n0) is A^4 + 4A^2 n0 + 2 n0^2
    m4 = float(np.mean(np.abs(seg) ** 4))
    rho = m4 / (s_pow ** 2 + 4 * s_pow * n0 + 2 * n0 ** 2)
    # spectrum: occupied bandwidth (99%) and power centroid
    nf = 1 << int(np.ceil(np.log2(max(len(seg), 256))))
    Sx = np.abs(np.fft.fftshift(np.fft.fft(seg, nf))) ** 2
    Sx = Sx - n0 * len(seg)
    Sx = np.clip(Sx, 0, None)
    fr = np.fft.fftshift(np.fft.fftfreq(nf, 1 / ofs))
    cs = np.cumsum(Sx) / max(Sx.sum(), 1e-30)
    f_l, f_h = fr[np.searchsorted(cs, 0.005)], fr[min(np.searchsorted(cs, 0.995), nf - 1)]
    centroid = float(np.sum(fr * Sx) / max(Sx.sum(), 1e-30))
    res = dict(t0_us=i0 / ofs * 1e6, t1_us=i1 / ofs * 1e6, f_center=fc + f_act + centroid,
               f_lo=fc + f_act + f_l, f_hi=fc + f_act + f_h, rho=rho, ofs=ofs, n0=n0,
               bw99=float(f_h - f_l))
    # in-band SNR (signal power over noise in the occupied band)
    res["snr_db"] = float(10 * np.log10(max(snr_slice * (2 * half) / max(f_h - f_l, 50e3), 1e-6)))
    res["const_env"] = bool(rho < 1.25) if snr_slice > 2 else None
    if snr_slice > 4 and len(seg) > 64:
        ifr = np.angle(seg[1:] * np.conj(seg[:-1])) * ofs / (2 * np.pi)
        R, conc, ph = _symbol_rate(ifr, ofs)
        res.update(sym_rate=R, sym_conc=conc)
        if R:
            sps = ofs / R
            v = ndimage.uniform_filter1d(ifr, max(1, int(sps / 2)))
            t0 = (ph / (2 * np.pi)) / R
            idx = ((t0 + 0.5 / R) * ofs + np.arange(int(len(v) / sps)) * sps).astype(int) % len(v)
            if len(idx) < 8:
                return res
            lo_, hi_ = _kmeans2(v[idx])
            res.update(dev_hz=float((hi_ - lo_) / 2), fsk_center=float(fc + f_act + (hi_ + lo_) / 2),
                       mod_index=float((hi_ - lo_) / R))
            if 0.85e6 < R < 1.15e6 and 2370e6 <= res["fsk_center"] <= 2520e6 and (i1 - i0) / ofs >= 120e-6:
                # BLE (2.4 GHz only -- widened by the +-20 MHz frequency-label errors that the
                # adverts themselves reveal -- and >= the shortest advert): re-slice +-0.75 MHz
                # at exactly 8 samples/bit around the FSK centre. Ungated, noise fragments with
                # a ~1 MHz "rate" made 20,000 such slices per 433 MHz block.
                yb, ob, _ = _slice(X, fs, res["fsk_center"] - fc, 0.75e6, 8e6)
                yw = yb[max(0, int(i0 / ofs * ob) - 80):int(i1 / ofs * ob) + 120]
                d = ble_try_decode(yw, ob, res["fsk_center"])
                if d:
                    res["ble"] = d
                else:
                    # XN297 toy drones / Bayang (research, PLAN R8): prefix 0x710F55 + Bayang checksum
                    d = XN.try_decode(yw, ob)
                    if d:
                        res["xn297"] = d
                    else:
                        # plain nRF24 toy drones (V2x2, SymaX; research, PLAN G6): fixed address + nRF CRC-16
                        d = TOYS.nrf24_try_decode(yw, ob)
                        if d:
                            res["xn297"] = d
                        else:
                            # Bluetooth BR/EDR (research, signal gaps S7): (64,30) BCH access code -> LAP
                            d = BR.try_decode(yw, ob)
                            if d:
                                res["bredr"] = d
            elif (0.45e6 < R < 0.55e6 and 2370e6 <= res["fsk_center"] <= 2520e6 and 0.6e6 <= res["bw99"] <= 2.0e6
                  and (i1 - i0) / ofs >= 400e-6):
                # BLE LE Coded S=8 (BT5 Long Range): the 0011 / 1100 patterns make a 1 Msym/s packet look like
                # 500 kBd at ~1 MHz width (research, PLAN T2b)
                yb, ob, _ = _slice(X, fs, res["fsk_center"] - fc, 0.75e6, 8e6)
                aa_bits = np.array([(BLE_ADV_AA >> i) & 1 for i in range(32)], np.uint8)
                d = ble_coded_decode(yb[max(0, int(i0 / ofs * ob) - 80):int(i1 / ofs * ob) + 120], ob,
                                     _ble_channel_index(res["fsk_center"]), aa_bits)
                if d:
                    res["ble"] = d
            elif 2370e6 <= res["fsk_center"] <= 2520e6 and 27e3 <= R <= 110e3 and 15e3 <= res["dev_hz"] <= 95e3:
                # FrSky D8 / D16 (CC2500; research, PLAN R6): a nominal rate within 8 % -> re-slice at 8 samples/bit,
                # sync 0xD391 / 0xD391D391 + length + CRC -> variant and TX ID
                if abs(R - RL.RATE) <= 0.08 * RL.RATE and (i1 - i0) / ofs * RL.RATE >= 32:
                    # RadioLink (CC2500; research, PLAN G4): doubled sync 0xABCD / 0xBADC / 0xABA5, PN9 dewhitening,
                    # CRC + checksum -> sub-protocol, TX ID, hop slot checked against the ID's hop table
                    yr, orr, _ = _slice(X, fs, res["fsk_center"] - fc, min(res["dev_hz"] + RL.RATE, 0.45 * 8 * RL.RATE),
                                        8 * RL.RATE)
                    j0 = max(0, int(i0 / ofs * orr) - 40)
                    j1 = max(int(i1 / ofs * orr) + 40, j0 + (8 + RL.MAX_LEN + 3) * 8 * 8 + 80)
                    d = RL.try_decode(yr[j0:j1], orr, res["fsk_center"])
                    if d:
                        res["radiolink"] = d
                kb = min(FRSKY.RATES.values(), key=lambda v: abs(v - R))
                if abs(R - kb) <= 0.08 * kb and (i1 - i0) / ofs * kb >= 64:
                    yr, orr, _ = _slice(X, fs, res["fsk_center"] - fc, min(res["dev_hz"] + kb, 0.45 * 8 * kb), 8 * kb)
                    # a weak packet splits into pieces: read from this piece's start over the longest packet
                    j0 = max(0, int(i0 / ofs * orr) - 40)
                    j1 = max(int(i1 / ofs * orr) + 40, j0 + (8 + 36) * 8 * 8 + 80)
                    d = FRSKY.try_identify(yr[j0:j1], orr, kb)
                    if d:
                        res["frsky"] = d
            elif res["fsk_center"] < 1.5e9 and res["dev_hz"] > 0:
                # SiK telemetry (research): an air rate within 5 % and deviation ~ rate (h ~ 2): re-slice at
                # 8 samples/bit around the FSK centre and look for the 0x2DD4 sync + MAVLink
                kb = min(SIK.RATES, key=lambda v: abs(v - R))
                if abs(R - kb) <= 0.05 * kb and 0.6 <= res["dev_hz"] / kb <= 1.4 and (i1 - i0) / ofs * kb >= 32:
                    half_s = min(res["dev_hz"] + kb, 0.45 * 8 * kb)
                    ys, os_, _ = _slice(X, fs, res["fsk_center"] - fc, half_s, 8 * kb)
                    # the energy detector can split a long packet at low SNR: decode from this piece's start over
                    # the longest SiK packet (preamble + sync + 5 + 255 bytes); pieces without the sync add nothing
                    j0 = max(0, int(i0 / ofs * os_) - 40)
                    j1 = max(int(i1 / ofs * os_) + 40, j0 + int((64 + 16 + 8 * 260) * 8 + 80))
                    d = SIK.sik_try_decode(ys[j0:j1], os_, kb)
                    if d:
                        d["crc_ok"] = bool(d.get("mavlink_crc_ok") or d.get("hw_crc_ok"))
                        res["sik"] = d
                if (abs(R - FL.CHIP_RATE) <= 0.06 * FL.CHIP_RATE and 30e3 <= res["dev_hz"] <= 70e3
                        and any(lo <= res["fsk_center"] <= hi for lo, hi in FL.BANDS)
                        and (i1 - i0) / ofs * FL.CHIP_RATE >= 40):
                    # >= 40 chips: a packet split into pieces is read from the piece holding preamble + sync
                    # FLARM beacon (research, PLAN T5): Manchester 100 kchip/s, sync F5 31 FA B6 -> ID, address type
                    yl, ol, _ = _slice(X, fs, res["fsk_center"] - fc, FL.DEV + FL.CHIP_RATE, 8 * FL.CHIP_RATE)
                    j0 = max(0, int(i0 / ofs * ol) - 40)
                    j1 = max(int(i1 / ofs * ol) + 40, j0 + (16 + 2 * 8 * 30) * 8 + 80)
                    # ADS-L SRD-860 M-band (research, PLAN G10): same PHY, sync F5 72 4B + length, CRC-24
                    da = ADSL.try_decode(yl[j0:max(j1, j0 + (16 + 2 * 8 * 32) * 8 + 80)], ol)
                    if da is None:
                        # OGN tracker (research, PLAN G10): same PHY, sync 0A F3 65 6C, Manchester-clean header
                        da = OGN.try_decode(yl[j0:max(j1, j0 + (16 + 2 * 8 * 30) * 8 + 80)], ol)
                    if da:
                        res["adsl"] = da
                    else:
                        d = FL.try_identify(yl[j0:j1], ol)
                        if d:
                            res["flarm"] = d
                if (abs(R - CF.RATE) <= 0.05 * CF.RATE and 30e3 <= res["dev_hz"] <= 55e3
                        and 850e6 <= res["fsk_center"] <= 935e6 and (i1 - i0) / ofs * CF.RATE >= 48):
                    # TBS Crossfire 150 Hz (research, PLAN R9): 0x3E6B3E sync after the preamble -> packet bytes
                    yc, oc, _ = _slice(X, fs, res["fsk_center"] - fc, CF.DEV + CF.RATE, 8 * CF.RATE)
                    # a weak packet splits into pieces: read from this piece's start over the longest packet
                    j0 = max(0, int(i0 / ofs * oc) - 40)
                    j1 = max(int(i1 / ofs * oc) + 40, j0 + (8 * 8 + 24 + 8 * 23) * 8 + 80)
                    d = CF.try_decode(yc[j0:j1], oc, res["fsk_center"])
                    if d:
                        res["crossfire"] = d
    R0 = res.get("sym_rate") or 0.0
    ble_branch = (0.85e6 < R0 < 1.15e6 and 2370e6 <= (res.get("fsk_center") or 0) <= 2520e6
                  and (i1 - i0) / ofs >= 120e-6)                        # XN297 was already tried there
    if (not (res.get("ble") or res.get("xn297") or ble_branch) and res.get("const_env") is not False
            and 2370e6 <= res["f_center"] <= 2520e6 and 100e-6 <= (i1 - i0) / ofs <= 600e-6
            and 0.3e6 <= res["bw99"] <= 1.6e6 and not 0.45e6 < R0 < 0.55e6):
        # XN297 toy drones at h ~0.3: the symbol-rate estimate is unreliable, so a packet-sized 2.4 GHz burst gets one
        # XN297 attempt at 1 Mb/s (prefix + Bayang checksum; research, PLAN R8)
        fx = res.get("fsk_center") or res["f_center"]
        yx, ox, _ = _slice(X, fs, fx - fc, 0.75e6, 8e6)
        yxs = yx[max(0, int(i0 / ofs * ox) - 80):int(i1 / ofs * ox) + 120]
        d = XN.try_decode(yxs, ox) or TOYS.nrf24_try_decode(yxs, ox)     # + plain nRF24 toys (research, PLAN G6)
        if d:
            res["xn297"] = d
    R = res.get("sym_rate")
    if (R and res.get("dev_hz") and not res.get("ble") and 0.45e6 < R < 0.55e6 and 2370e6 <= res["fsk_center"] <= 2520e6
            and 120e3 <= res["dev_hz"] <= 260e3 and (i1 - i0) / ofs * R >= 64):
        # FlySky AFHDS 2A (A7105; research, PLAN R7): ID code 0x5475C52A -> identification; FEC + CRC under candidate
        # conventions -> TX / RX IDs. Read from this piece's start over the longest packet
        ya, oa, _ = _slice(X, fs, res["fsk_center"] - fc, AF2.DEV + AF2.RATE, 8 * AF2.RATE)
        j0 = max(0, int(i0 / ofs * oa) - 40)
        j1 = max(int(i1 / ofs * oa) + 40, j0 + (64 + (AF2.PAYLOAD + 2) * 14) * 8 + 80)
        d = AF2.try_identify(ya[j0:j1], oa, f_hz=res["fsk_center"])
        if d:
            res["afhds2a"] = d
    if (R and res.get("dev_hz") and not res.get("ble") and 250e3 <= R <= 350e3 and 60e3 <= res["dev_hz"] <= 140e3
            and (2400e6 <= res["fsk_center"] <= 2480e6 or 860e6 <= res["fsk_center"] <= 930e6)
            and (i1 - i0) / ofs <= 700e-6):
        # ExpressLRS LR1121 GFSK modes (research, PLAN G3): a SYNC packet whose CRC verifies with its own UID bytes,
        # which must equal the 16-bit sync word -> UID4 / UID5, OTA version, whitening variant
        ye, oe, _ = _slice(X, fs, res["fsk_center"] - fc, ELRS.GFSK_DEV + ELRS.GFSK_RATE, 8 * ELRS.GFSK_RATE)
        d = ELRS.gfsk_try_decode(ye[max(0, int(i0 / ofs * oe) - 80):int(i1 / ofs * oe) + 120], oe)
        if d:
            res["elrs_gfsk"] = d
    if (R and res.get("dev_hz") and not res.get("ble") and 2400e6 <= res["fsk_center"] <= 2490e6
            and 120e3 <= res["dev_hz"] <= 260e3 and 1.5e-3 <= (i1 - i0) / ofs <= 2.6e-3):
        # Hubsan (A7105 100 kb/s; research, PLAN G5): session ID code after the preamble, checksum + model bytes.
        # At h ~3.7 the symbol-rate estimate is unreliable (35-100 kHz): gated on deviation and packet length
        # (1.9-2.1 ms) instead; the decode verifies itself
        yh, oh, _ = _slice(X, fs, res["fsk_center"] - fc, HB.DEV + HB.RATE, 8 * HB.RATE)
        j0 = max(0, int(i0 / ofs * oh) - 40)
        j1 = max(int(i1 / ofs * oh) + 40, j0 + (64 + 8 * (HB.LEN + 2)) * 8 + 80)
        d = HB.try_decode(yh[j0:j1], oh, res["fsk_center"])
        if d:
            res["hubsan"] = d
    dur = (i1 - i0) / ofs
    if (R and res.get("dev_hz") and not res.get("ble") and 2400e6 <= res["fsk_center"] <= 2490e6
            and 200e3 <= R <= 310e3 and 30e3 <= res["dev_hz"] <= 110e3 and 0.6e-3 <= dur <= 4.4e-3):
        # Graupner HoTT (MSK + CC2500 FEC) and Hitec Optima / Minima (MSK, sync = TX ID; research, PLAN G12):
        # 249.94 kBd, h = 0.5; a weak packet splits into pieces, so each is read over the longest packet
        yk, ok_, _ = _slice(X, fs, res["fsk_center"] - fc, HOTT.DEV + HOTT.RATE, 8 * HOTT.RATE)
        j0 = max(0, int(i0 / ofs * ok_) - 80)
        d = HOTT.try_decode(yk[j0:max(int(i1 / ofs * ok_) + 80, j0 + (64 + HOTT.coded_len(HOTT.LEN_TX)) * 8 + 160)],
                            ok_, res["fsk_center"])
        if d is None:
            d = HITEC.hitec_try_decode(yk[j0:max(int(i1 / ofs * ok_) + 80, j0 + (8 * 60) * 8 + 160)], ok_,
                                       res["fsk_center"])
        if d:
            res["cc2500_rc"] = d
    if (R and res.get("dev_hz") and not res.get("ble") and 2400e6 <= res["fsk_center"] <= 2450e6
            and 100e3 <= R <= 160e3 and 20e3 <= res["dev_hz"] <= 60e3 and 0.4e-3 <= dur <= 1.9e-3):
        # Futaba S-FHSS (CC2500 2-FSK 128 kBd, sync D391, CRC in CC2400 mode; research, PLAN G12)
        ys, os_, _ = _slice(X, fs, res["fsk_center"] - fc, HITEC.SF_DEV + HITEC.SF_RATE, 8 * HITEC.SF_RATE)
        j0 = max(0, int(i0 / ofs * os_) - 80)
        d = HITEC.sfhss_try_decode(ys[j0:max(int(i1 / ofs * os_) + 80, j0 + 8 * 8 * 24 + 160)], os_, res["fsk_center"])
        if d:
            res["cc2500_rc"] = d
    if (R and res.get("dev_hz") and not (res.get("ble") or res.get("sik") or res.get("crossfire")
                                          or res.get("frsky") or res.get("afhds2a") or res.get("xn297")
                                          or res.get("radiolink") or res.get("hubsan") or res.get("elrs_gfsk")
                                          or res.get("cc2500_rc") or res.get("flarm") or res.get("adsl"))
            and 20e3 <= R <= 1.2e6 and res["snr_db"] >= 10 and (i1 - i0) / ofs * R >= FF.MIN_PRE + FF.SYNC_BITS + 16):
        # any other FSK packet: blind frame structure (preamble, 32-bit sync-word candidate, length; research)
        half_f = min(res["dev_hz"] + R, 0.45 * 8 * R)
        yf, of_, _ = _slice(X, fs, res["fsk_center"] - fc, half_f, 8 * R)
        fr = FF.frame(yf[max(0, int(i0 / ofs * of_) - 16):int(i1 / ofs * of_) + 16], of_, R)
        if fr:
            res["frame"] = fr
    if (not res.get("sym_rate") and 2370e6 <= res["f_center"] <= 2520e6 and (i1 - i0) / ofs >= 120e-6
            and 0.5e6 <= res["bw99"] <= 2.5e6):
        # low SNR (PLAN L1): no symbol rate below ~6 dB slice SNR, but a BLE-wide, advert-long 2.4 GHz burst is
        # worth one coherent access-address search + decode (gfsk_coherent; nothing is reported without the CRC)
        yb, ob, _ = _slice(X, fs, res["f_center"] - fc, 0.75e6, 8e6)
        aa_bits = np.array([(BLE_ADV_AA >> i) & 1 for i in range(32)], np.uint8)
        ys = yb[max(0, int(i0 / ofs * ob) - 80):int(i1 / ofs * ob) + 120]
        chx = _ble_channel_index(res["f_center"])
        d = ble_coherent_decode(ys, ob, chx, aa_bits) or ble_coded_decode(ys, ob, chx, aa_bits)
        if d:
            res["ble"] = d
    return res


def _gsm_burst(r, f):
    """One to eight GSM time slots on the 200 kHz raster in a GSM band: a GSM carrier, not an RC link."""
    dur = r["t1_us"] - r["t0_us"]
    k = round(dur / GSM_SLOT_US)
    if not 1 <= k <= 8 or abs(dur - k * GSM_SLOT_US) > GSM_SLOT_TOL_US:
        return False
    if not any("gsm_ul" in b[4] or "gsm_dl" in b[4] for b in B.lookup(f)):
        return False
    return abs(f - 200e3 * round(f / 200e3)) <= GSM_RASTER_TOL_HZ


def _fsk_classify(r, fc=2.44e9):
    """(family, label, confidence, protocol nominal symbol rate or None)."""
    if r.get("esb") and not r.get("ble"):
        return "esb_link", "esb_like_2m", 0.6, None
    if r.get("ieee802154"):
        ok = r["ieee802154"]["crc_ok"]
        return "telemetry", "ieee802154" if ok else "ieee802154_like", 0.9 if ok else 0.6, None
    if r.get("const_env") is False:
        return ("unknown", "nb_nonconst_envelope", 0.4, None) \
            if r.get("snr_db", 0) >= 12 and r["t1_us"] - r["t0_us"] >= 60 else (None, None, 0, None)
    R, dev = r.get("sym_rate"), r.get("dev_hz")
    if (r.get("sik") or {}).get("crc_ok"):
        return "telemetry", "sik_telemetry", 0.9, r["sik"]["air_rate_bps"]
    if r.get("crossfire"):
        return "rc_fhss", "crossfire_sync", 0.8, CF.RATE
    if r.get("adsl"):
        if r["adsl"]["protocol"].startswith("OGN"):
            return "telemetry", "ogn_like", 0.75, FL.CHIP_RATE
        return "telemetry", "adsl", 0.9, FL.CHIP_RATE
    if r.get("flarm"):
        return "telemetry", "flarm_like", 0.8, FL.CHIP_RATE
    if r.get("xn297"):
        d = r["xn297"]
        if d.get("payload_len"):                          # a named toy protocol (toys.py)
            return "rc_fhss", "toy_drone", 0.85 if d.get("checksum_ok") else 0.6, XN.RATE
        ok = d.get("checksum_ok")
        return "rc_fhss", "bayang" if ok else "xn297_like", 0.85 if ok else 0.6, XN.RATE
    if r.get("elrs_gfsk"):
        return "rc_fhss", "elrs_gfsk", 0.9, ELRS.GFSK_RATE
    if r.get("radiolink"):
        return "rc_fhss", "radiolink", 0.9, RL.RATE
    if r.get("hubsan"):
        return "rc_fhss", "hubsan", 0.85, HB.RATE
    if r.get("cc2500_rc"):
        p = r["cc2500_rc"]["protocol"]
        if p.startswith("Graupner"):
            return "rc_fhss", "hott", 0.9, HOTT.RATE
        if p.startswith("Hitec"):
            return "rc_fhss", "hitec", 0.9, HITEC.RATE
        return "rc_fhss", "sfhss", 0.85, HITEC.SF_RATE
    if r.get("afhds2a"):
        if r["afhds2a"]["protocol"].startswith("FlySky AFHDS 1"):
            return "rc_fhss", "afhds1" if r["afhds2a"]["crc_ok"] else "afhds1_like", \
                0.85 if r["afhds2a"]["crc_ok"] else 0.7, AF2.RATE
        return "rc_fhss", "afhds2a" if r["afhds2a"]["crc_ok"] else "afhds2a_id", \
            0.85 if r["afhds2a"]["crc_ok"] else 0.7, AF2.RATE
    if r.get("frsky"):
        d = r["frsky"]
        nom = FRSKY.RATES.get(d["variant"].replace(" LBT v1", " LBT"), FRSKY.RATES["D16 v1"])
        return "rc_fhss", "frsky_" + d["variant"].split()[0].lower(), 0.9 if d["crc_ok"] else 0.7, nom
    if r.get("ble"):
        return "bluetooth", "ble_adv" if r["ble"].get("crc_ok") else "bluetooth_like", \
            0.95 if r["ble"].get("crc_ok") else 0.8, 1e6
    if r.get("bredr"):                                    # BCH-verified access code (research, signal gaps S7)
        return "bluetooth", "bluetooth_bredr", 0.9, 1e6
    # protocol-level labels need a clean measurement: constant envelope, a
    # clear symbol clock (transition-phase concentration) and >= 10 dB
    if not R or not dev or r.get("sym_conc", 0) < 0.5 or r.get("snr_db", 0) < 6 or not r.get("const_env"):
        return "narrowband", "nb_burst", 0.3, None
    # and room for a packet: preamble + sync + header alone are >= 64 symbols in every
    # family below (a 90 us blip at 180 kBd is 16 symbols -- not a SiK/Crossfire packet).
    # Fragments of a weak packet (9-22 symbols each) cannot be told from such blips by
    # rate agreement or timing (tried on the 433 MHz artifacts), so they stay nb_burst:
    # FrSky-like is then labelled from 12.5 dB instead of 7.5 dB, detected from 5 dB.
    if (r["t1_us"] - r["t0_us"]) * 1e-6 * R < 64:
        return "narrowband", "nb_burst", 0.3, None
    if 0.85e6 < R < 1.15e6 and 100e3 < dev < 400e3:
        return "bluetooth", "bluetooth_like", 0.6, 1e6      # BLE 1M (h~0.5) or BR (h~0.32) or DSMX-like
    carson = 2 * (dev + R / 2)
    if r["bw99"] > 3 * carson + 100e3:                    # spectrum far wider than the FSK explains
        return "narrowband", "fsk_unknown", 0.3, None
    sub_ghz = {"crossfire_fsk_like", "sik_like"}
    for lab, (r0, r1), (d0, d1) in FSK_CLASSES:
        if (lab in sub_ghz) != (fc < 1.5e9):
            continue
        if r0 <= R <= r1 and d0 <= dev <= d1:
            nom = {"elrs_flrc_like": 650e3, "elrs_gfsk_like": 300e3, "afhds2a_like": 500e3,
                   "crossfire_fsk_like": 85.1e3, "frsky_d16_like": 70.0e3, "sik_like": None}[lab]
            return "rc_fhss", lab, 0.6, nom
    return "narrowband", "fsk_unknown", 0.4, None


def _local_fft(x, fs, t0_us, t1_us, margin_us):
    """FFT of x around [t0, t1] (+-margin) and the segment start in us. A burst needs
    only its own milliseconds: slicing the whole capture's FFT per burst cost minutes
    when artifacts produced thousands of candidates."""
    i0 = max(0, int((t0_us - margin_us) * 1e-6 * fs))
    i1 = min(len(x), int((t1_us + margin_us) * 1e-6 * fs) + 1)
    return sfft.fft(np.asarray(x[i0:i1], np.complex64), workers=-1), i0 / fs * 1e6


def _esb_try(x, fs, fc, r):
    """2 Mb/s ESB-like packet in a measured 2.4 GHz burst 0.5-3 MHz wide: slice +-1.25 MHz at 8 samples/bit from
    an FFT of the burst +-50 us only (the slice's filter edges settle), cropped to the burst -2 / +6 us (neighbouring
    packets would bias the bit threshold and timing); kept with a >= 16-bit 0101 preamble (esb.slice_packet)."""
    if not (ESB_BAND[0] <= r["f_center"] <= ESB_BAND[1] and ESB_US[0] <= r["t1_us"] - r["t0_us"] <= ESB_US[1]
            and ESB_BW[0] <= r["bw99"] <= ESB_BW[1] and r["snr_db"] >= 6.0):
        return None
    Xs, off = _local_fft(x, fs, r["t0_us"], r["t1_us"], ESB_MARGIN_US)
    y, oe, fe = _slice(Xs, fs, r["f_center"] - fc, ESB.HALF_BW, ESB.OUT_FS)
    i0 = max(0, int((r["t0_us"] - off - 2.0) * 1e-6 * oe))
    d = ESB.slice_packet(y[i0:int((r["t1_us"] - off + 6.0) * 1e-6 * oe)], oe)
    if d:
        d["f_carrier"] = fc + fe + d["f_off_hz"]
    return d


def _zb_try(x, fs, fc, r):
    """IEEE 802.15.4 O-QPSK frame (research, PLAN T4) in a 2.4 GHz burst on an 802.15.4 channel, 1.2-4 MHz wide and
    >= 250 us (SHR + PHR + a short MPDU): slice +-1.6 MHz around the nominal channel at 8 MHz (4 samples / chip)
    from a local FFT of the burst +-50 us, cropped to the burst +-20 us."""
    ch = ZB.channel(r["f_center"], tol_hz=300e3)
    if ch is None or not 1.2e6 <= r["bw99"] <= 4e6 or r["t1_us"] - r["t0_us"] < 250 or r["snr_db"] < 3.0:
        return None
    f0 = 2405e6 + 5e6 * (ch - 11)
    Xs, off = _local_fft(x, fs, r["t0_us"], r["t1_us"], ESB_MARGIN_US)
    y, oz, _ = _slice(Xs, fs, f0 - fc, 1.6e6, 8e6)
    i0 = max(0, int((r["t0_us"] - off - 20.0) * 1e-6 * oz))
    d = ZB.demod(y[i0:int((r["t1_us"] - off + 20.0) * 1e-6 * oz)], oz)
    if d:
        d["channel"] = ch
    return d


_CARRY_KEYS = ("xn297", "radiolink", "hubsan", "elrs_gfsk", "afhds2a", "cc2500_rc")


def _refine_local(x, fs, fc, cand, noise_win_us=2000.0):
    """fsk_refine on a segment around the candidate; times come back capture-relative."""
    Xs, off = _local_fft(x, fs, cand["t0_us"], cand["t1_us"], noise_win_us + 500.0)
    rs = fsk_refine(Xs, fs, fc, dict(cand, t0_us=cand["t0_us"] - off, t1_us=cand["t1_us"] - off), noise_win_us)
    for r in rs:
        r["t0_us"] += off; r["t1_us"] += off
    return rs


def detect_fsk(x, X, fs, fc, log=None, **kw):
    cands = fsk_scan(x, fs, fc, **kw)
    if log: log(f"  fsk: {len(cands)} narrowband candidates")
    refined = [r for c in cands for r in _refine_local(x, fs, fc, c)]
    # fragments of one long low-SNR burst: merge same-frequency pieces separated
    # by < 100 us and re-measure the union
    refined.sort(key=lambda r: r["t0_us"])
    merged = []
    for r in refined:
        for m in merged:
            if abs(m["f_center"] - r["f_center"]) < 0.5 * max(m["bw99"], r["bw99"], 100e3) and \
                    r["t0_us"] - m["t1_us"] < 100 and r["t1_us"] - m["t0_us"] > -100:
                m["parts"].append(r); m["t1_us"] = max(m["t1_us"], r["t1_us"]); m["t0_us"] = min(m["t0_us"], r["t0_us"])
                break
        else:
            merged.append(dict(r, parts=[r]))
    final = []
    for m in merged:
        if len(m["parts"]) == 1:
            final.append(m["parts"][0]); continue
        lo = min(p["f_lo"] for p in m["parts"]); hi = max(p["f_hi"] for p in m["parts"])
        rr = _refine_local(x, fs, fc, dict(t0_us=m["t0_us"], t1_us=m["t1_us"], f_lo=lo, f_hi=hi, fc_b=(lo + hi) / 2))
        # a decode found on one piece (the new PLAN G3-G6 decoders) stays with the merged burst that covers it
        for part in m["parts"]:
            for k in _CARRY_KEYS:
                if part.get(k) and (k != "xn297" or part[k].get("payload_len")) and \
                        (k != "afhds2a" or part[k]["protocol"].startswith("FlySky AFHDS 1")):
                    host = next((q for q in rr if q["t0_us"] - 50 <= part["t0_us"] <= q["t1_us"]
                                 and not any(q.get(kk) for kk in _CARRY_KEYS)), None)
                    if host is not None:
                        host[k] = part[k]
        final += rr if rr else m["parts"]
    out, esb, covered = [], [], []
    for r in final:
        # a decoded long packet (RadioLink, 9.3 ms) splits into pieces: later pieces inside its air time join it
        f0 = r.get("fsk_center", r["f_center"])
        host = next((c for c in covered if abs(c[0] - f0) < 100e3 and c[1] <= r["t0_us"] < c[2]), None)
        if host is not None:
            out[host[3]].end_us = max(out[host[3]].end_us, r["t1_us"])
            continue
        # energy significance of the refined burst: snr * sqrt(B*T) (independent
        # time-frequency cells); noise blobs that survive the scan sit at 1-3
        sig = 10 ** (r["snr_db"] / 10) * np.sqrt(max(r["bw99"], 50e3) * (r["t1_us"] - r["t0_us"]) * 1e-6)
        if r["t1_us"] - r["t0_us"] < 20 or sig < 8.0:
            continue
        f = r.get("fsk_center", r["f_center"])
        if _gsm_burst(r, f):
            continue
        z = _zb_try(x, fs, fc, r)
        if z:
            r["ieee802154"] = z
        d = None if z else _esb_try(x, fs, fc, r)
        if d:
            r["esb"] = d
        fam, lab, conf, nom = _fsk_classify(r, fc)
        if fam is None:                                   # weak non-constant-envelope piece: not ours
            continue
        feats = dict(modulation="gfsk" if r.get("dev_hz") else None, bw_hz=r["bw99"], bw_mhz=r["bw99"] / 1e6,
                     burst_us=r["t1_us"] - r["t0_us"], symbol_rate=r.get("sym_rate"), dev_hz=r.get("dev_hz"),
                     mod_index=r.get("mod_index"), const_env_rho=r["rho"], sym_conc=r.get("sym_conc"))
        if nom and r.get("sym_rate"):
            feats["clock_ppm"] = (r["sym_rate"] / nom - 1) * 1e6
            feats["symbol_us"] = 1e6 / r["sym_rate"]
        if fam == "bluetooth" and r.get("ble"):
            fn = 2402e6 + 2e6 * round((f - 2402e6) / 2e6)
            feats["cfo_hz"] = f - fn; feats["cfo_ppm"] = (f - fn) / fn * 1e6
        if fam == "esb_link":
            d = r["esb"]; f = d["f_carrier"]; fn = 1e6 * round(f / 1e6)       # ESB channels: 2400 + k MHz
            feats.update(modulation="gfsk", symbol_rate=ESB.RATE, symbol_us=1e6 / ESB.RATE, dev_hz=d["dev_hz"],
                         mod_index=2 * d["dev_hz"] / ESB.RATE, sym_conc=None, preamble_bits=d["preamble_bits"],
                         cfo_hz=f - fn, cfo_ppm=(f - fn) / fn * 1e6)
            esb.append((len(out), d["bits"]))
        if r.get("frame") and fam != "esb_link":
            feats["frame"] = r["frame"]
        out.append(Emission(family=fam, label=lab, start_us=r["t0_us"], end_us=r["t1_us"],
                            f_lo_mhz=min(r["f_lo"], f - r["bw99"] / 2) / 1e6,
                            f_hi_mhz=max(r["f_hi"], f + r["bw99"] / 2) / 1e6, center_mhz=f / 1e6,
                            snr_db=r["snr_db"], confidence=conf, decoded=r.get("ble") or r.get("bredr") or r.get("sik") or r.get("crossfire") or r.get("frsky")
                            or r.get("afhds2a") or r.get("xn297") or r.get("ieee802154") or r.get("flarm") or r.get("adsl")
                            or r.get("radiolink") or r.get("hubsan") or r.get("elrs_gfsk") or r.get("cc2500_rc"),
                            features=feats))
        if r.get("radiolink"):
            n_bits = 8 * (4 + 4 + r["radiolink"]["length"] + 3)
            covered.append((f, r["t0_us"], r["t0_us"] + n_bits / RL.RATE * 1e6 + 200.0, len(out) - 1))
        if r.get("adsl"):                                 # ADS-L / OGN (FLARM PHY): the packet's chips after the sync
            n_chips = 2 * 8 * ((r["adsl"].get("length") or 26) + 5) if not r["adsl"]["protocol"].startswith("OGN") \
                else 2 * 8 * (4 + 26)
            covered.append((f, r["t0_us"], r["t0_us"] + n_chips / FL.CHIP_RATE * 1e6 + 200.0, len(out) - 1))
    if esb:                                       # one CRC search per link (packets sharing a sync word)
        for (i, _), v in zip(esb, ESB.verify_links([b for _, b in esb], [out[i].start_us for i, _ in esb])):
            out[i].decoded = v
            if v["crc_ok"]:
                out[i].confidence = 0.85
    return out


# ============================================================== combine
def suppress_images(ems, fc, min_db=3.0):
    """Drop bursts that are the IQ-imbalance mirror (2*fc - f) of a stronger,
    time-overlapping burst (Pluto/AD936x image rejection is only ~10-40 dB)."""
    keep = []
    cm = np.array([o.center_mhz for o in ems])
    for i, e in enumerate(ems):
        fm = 2 * fc / 1e6 - e.center_mhz
        bw = max(e.f_hi_mhz - e.f_lo_mhz, 0.2)
        img = False
        for j in np.flatnonzero(np.abs(cm - fm) <= bw / 2 + 0.1):   # only bursts near the mirror frequency
            o = ems[j]
            if j == i:
                continue
            ov = min(o.end_us, e.end_us) - max(o.start_us, e.start_us)
            if ov > 0.5 * (e.end_us - e.start_us) and (o.snr_db or -99) > (e.snr_db or -99) + min_db and \
                    (e.end_us - e.start_us) <= 1.3 * (o.end_us - o.start_us):
                img = True; break
        if not img:
            keep.append(e)
    return keep


def detect(x, fs, fc, ctx=None):
    """All RC / telemetry bursts in x (complex baseband at fs, centre fc).
    ctx keys used: 'rc_fft' (cached FFT of x), 'log' (callable), 'lora_bws', 'lora_only'
    ((lo, hi) Hz offsets from fc: only the LoRa scan, on sub-bands centred there; no FSK/BLE/ESB)."""
    ctx = ctx if ctx is not None else {}
    log = ctx.get("log")
    X = ctx.get("rc_fft")
    if X is None or len(X) != len(x):
        X = sfft.fft(np.asarray(x, np.complex64), workers=-1)
        ctx["rc_fft"] = X
    ems = detect_lora(X, fs, fc, bws=ctx.get("lora_bws"), log=log, f_range=ctx.get("lora_only"))
    fsk = [] if ctx.get("lora_only") else detect_fsk(x, X, fs, fc, log=log)
    # an FSK candidate covering a LoRa packet is the same burst seen by energy
    fsk = [e for e in fsk if not any(abs(e.center_mhz - l.center_mhz) < 0.6 * (l.f_hi_mhz - l.f_lo_mhz) and
                                     min(e.end_us, l.end_us) - max(e.start_us, l.start_us) > 0
                                     for l in ems)]
    ems = suppress_images(ems + fsk, fc)
    ems = [e for e in ems if not _is_raw_image(x, fs, fc, e)]
    ems = sorted(ems, key=lambda e: e.start_us)
    _elrs_hops(ems)
    return ems


def _elrs_hops(ems):
    """ExpressLRS hop recovery (research, PLAN R2d): a CRC-verified SYNC packet plus the link's other packets in the
    block (same SF and bandwidth) -> UID2 / UID3 candidates and, when unique, the next channels; attached to the SYNC
    packet's decoded dict as hop_recovery. Overlapping domain grids (IN866 inside EU868, AU915 inside FCC915) are all
    tried with the packets on their grid (off-grid emissions are not the link); a domain is kept when the search leaves
    candidates."""
    for s in ems:
        d = s.decoded or {}
        if d.get("protocol") != "ExpressLRS" or d.get("packet_type") != "SYNC" or not d.get("crc_ok"):
            continue
        modes = s.features.get("elrs_modes") or []
        if not modes:
            continue
        _, _, _, _, _, _, hop, ival = ELRS_LORA[modes[0]]
        link = [e for e in ems if e is not s and e.features.get("modulation") == "lora"
                and e.features.get("sf") == s.features.get("sf") and e.features.get("bw_hz") == s.features.get("bw_hz")]
        fits = []
        for dom, (f0, f1, n) in ELRS_DOMAINS.items():
            sp = (f1 - f0) / (n - 1)
            idx = lambda f: int(round((f - f0) / sp))                                   # noqa: E731
            ok = lambda f: 0 <= idx(f) < n and abs(f - (f0 + idx(f) * sp)) < 150e3     # noqa: E731
            fs_ = s.features.get("f_carrier_hz")
            if fs_ is None or not ok(fs_):
                continue
            pk = [(e.start_us, idx(e.features["f_carrier_hz"])) for e in link if ok(e.features["f_carrier_hz"])]
            r = ELRS.recover_hops(d, s.start_us, idx(fs_), ival, hop, n, pk)
            if r and r["uid23_candidates"]:
                fits.append(dict(r, mode=modes[0], domain=dom))
        if fits:
            d["hop_recovery"] = dict(fits[0], other_domains=[f["domain"] for f in fits[1:]])


def _is_raw_image(x, fs, fc, e, min_db=3.0):
    """Mirror test on the raw signal (catches images whose parent burst was
    not itself detected, e.g. hidden under a wideband video burst)."""
    bw = max((e.f_hi_mhz - e.f_lo_mhz) * 1e6, 200e3)
    ofs = float(min(max(2 * bw, 1e6), 8e6, fs))
    X, off = _local_fft(x, fs, e.start_us, e.end_us, 50.0)
    a, b = (e.start_us - off) * 1e-6, (e.end_us - off) * 1e-6
    p = []
    for f in (e.center_mhz * 1e6, 2 * fc - e.center_mhz * 1e6):
        y, o, _ = _slice(X, fs, f - fc, bw / 2, ofs)
        p.append(np.mean(np.abs(y[int(a * o):max(int(b * o), int(a * o) + 1)]) ** 2))
    return 10 * np.log10(p[1] / max(p[0], 1e-30)) > min_db


# ============================================================== hop-set grouping
def _theil_sen(n, t):
    """Median of pairwise slopes: one burst with a one-symbol timing slip
    does not bend the packet period."""
    sl = [(t[j] - t[i]) / (n[j] - n[i]) for i in range(len(n)) for j in range(i + 1, len(n)) if n[j] != n[i]]
    return float(np.median(sl))


def group_hops(emissions, t_tol_us=None, ppm_tol=3.0, min_members=3, max_period_us=50000.0):
    """RC-specific link grouping. Bursts with identical PHY parameters (LoRa
    SF/BW/IQ, or FSK class + symbol rate) whose start times fit one packet grid
    t = t0 + n*T (missed packets allowed) form one link; cfo_ppm (one crystal
    -> one ppm across hops) splits links sharing a grid. Returns list of dict:
    members (indices), packet_interval_us, packet_rate_hz, timer_ppm (vs the
    nearest known RC interval), channels_mhz (in time order), cfo_ppm (median),
    label, protocol_guess."""
    def key(e):
        f = e.features
        if f.get("modulation") == "lora":
            return ("lora", f["sf"], round(f["bw_hz"]), bool(f.get("iq_inverted")))
        if e.family in ("rc_fhss", "narrowband", "bluetooth") and f.get("symbol_rate"):
            return ("fsk", e.label, round(f["symbol_rate"] / 20e3))
        return None
    pools = {}
    for i, e in enumerate(emissions):
        k = key(e)
        if k is not None:
            pools.setdefault(k, []).append(i)
    links = []
    known = sorted({v for vs in KNOWN_INTERVALS_US.values() for v in vs})
    for k, idx in pools.items():
        remaining = sorted(idx, key=lambda i: emissions[i].start_us)
        while len(remaining) >= min_members:
            t = np.array([emissions[i].start_us for i in remaining])
            ppm = np.array([emissions[i].features.get("cfo_ppm") if emissions[i].features.get("cfo_ppm") is not None
                            else np.nan for i in remaining])
            tol = t_tol_us or max(5.0, 0.02 * np.median([emissions[i].end_us - emissions[i].start_us for i in remaining]))
            best = None
            # candidate periods: all pairwise differences / small integers, plus known RC intervals
            cand_T = set(known)
            for a in range(len(t)):
                for b in range(a + 1, min(len(t), a + 6)):
                    for m in (1, 2, 3, 4):
                        d = (t[b] - t[a]) / m
                        if 400 <= d <= max_period_us:
                            cand_T.add(round(d, 1))
            for T in sorted(cand_T):
                for a in range(len(t)):
                    n = np.round((t - t[a]) / T)
                    res = t - t[a] - n * T
                    m = np.abs(res) <= tol
                    if np.isfinite(ppm[a]):
                        m &= ~(np.abs(ppm - ppm[a]) > ppm_tol)       # nan ppm -> allowed
                    if m.sum() < min_members:
                        continue
                    # least-squares period over members; prefer more members, then larger T
                    nn, tt = n[m], t[m]
                    if np.ptp(nn) == 0:
                        continue
                    T_fit = np.polyfit(nn, tt, 1)[0]
                    # fill = occupied fraction of grid slots over the span: a
                    # sub-multiple period (T/2 interleaving two same-rate links)
                    # fits as many bursts but leaves half the slots empty
                    fill = m.sum() / (np.ptp(nn) + 1)
                    score = (round(m.sum() * min(fill, 1.0), 3), T)
                    if best is None or score > best[0]:
                        best = (score, m.copy(), T_fit, a)
            if best is None:
                break
            _, m, T_fit, a = best
            if np.isfinite(ppm[m]).sum() >= 2:                # drop members off the link's crystal ppm
                med = np.nanmedian(ppm[m])
                m &= ~(np.abs(ppm - med) > ppm_tol)
            mem = [remaining[i] for i in np.flatnonzero(m)]
            # period from SFD-verified members only when possible (preamble-only
            # detections carry +-1 symbol timing)
            nn = np.round((t[m] - t[a]) / T_fit)
            good = np.isfinite(ppm[m])
            if good.sum() >= 3 and np.ptp(nn[good]) > 0:
                T_fit = _theil_sen(nn[good], t[m][good])
            Tn = min(known, key=lambda v: abs(v - T_fit))
            if good.sum() < 3 and abs(T_fit / Tn - 1) < 5e-3:
                T_fit = float(Tn)                         # too few precise starts: snap to the known RC interval
            e0 = emissions[mem[0]]
            ppms = [emissions[i].features.get("cfo_ppm") for i in mem if emissions[i].features.get("cfo_ppm") is not None]
            guess = None
            if k[0] == "lora":
                modes = [mm for mm, v in ELRS_LORA.items() if v[1] == k[1] and abs(v[2] - k[2]) < 1e3 and abs(v[7] - T_fit) < 0.01 * v[7]]
                guess = ("expresslrs_" + modes[0]) if modes else None
            links.append(dict(members=mem, label=e0.label, family=e0.family, packet_interval_us=float(T_fit),
                              packet_rate_hz=float(1e6 / T_fit),
                              timer_ppm=float((T_fit / Tn - 1) * 1e6) if abs(T_fit / Tn - 1) < 300e-6 else None,
                              channels_mhz=[round(emissions[i].center_mhz, 3) for i in mem],
                              n_channels=len({round(emissions[i].center_mhz, 1) for i in mem}),
                              cfo_ppm=float(np.median(ppms)) if ppms else None,
                              cfo_ppm_spread=float(np.std(ppms)) if len(ppms) > 1 else None,
                              protocol_guess=guess))
            remaining = [i for i in remaining if i not in mem]
    # merge pieces of one link split by missed packets: same PHY key/label, same
    # period (0.2%), same crystal ppm, start times on one common grid
    merged = True
    while merged:
        merged = False
        for i in range(len(links)):
            for j in range(i + 1, len(links)):
                A, B = links[i], links[j]
                if A["label"] != B["label"] or abs(A["packet_interval_us"] / B["packet_interval_us"] - 1) > 2e-3:
                    continue
                if A["cfo_ppm"] is not None and B["cfo_ppm"] is not None and abs(A["cfo_ppm"] - B["cfo_ppm"]) > ppm_tol:
                    continue
                T = A["packet_interval_us"]
                d = emissions[B["members"][0]].start_us - emissions[A["members"][0]].start_us
                if abs(d - round(d / T) * T) > max(10.0, 0.05 * T):
                    continue
                mem = sorted(A["members"] + B["members"], key=lambda k: emissions[k].start_us)
                ts = np.array([emissions[k].start_us for k in mem]); nn = np.round((ts - ts[0]) / T)
                Tf = _theil_sen(nn, ts) if np.ptp(nn) > 0 else T
                ppms = [emissions[k].features.get("cfo_ppm") for k in mem if emissions[k].features.get("cfo_ppm") is not None]
                Tn = min(known, key=lambda v: abs(v - Tf))
                A.update(members=mem, packet_interval_us=Tf, packet_rate_hz=1e6 / Tf,
                         timer_ppm=float((Tf / Tn - 1) * 1e6) if abs(Tf / Tn - 1) < 300e-6 else None,
                         channels_mhz=[round(emissions[k].center_mhz, 3) for k in mem],
                         n_channels=len({round(emissions[k].center_mhz, 1) for k in mem}),
                         cfo_ppm=float(np.median(ppms)) if ppms else None,
                         cfo_ppm_spread=float(np.std(ppms)) if len(ppms) > 1 else None,
                         protocol_guess=A["protocol_guess"] or B["protocol_guess"])
                links.pop(j); merged = True
                break
            if merged:
                break
    return links
