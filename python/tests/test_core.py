"""Fast unit tests for the pieces with exact, checkable behaviour."""
import numpy as np

from spectrum_intel.receivers.adsb import modes as adsb
from spectrum_intel.receivers.gsm import sch as gsm
from spectrum_intel.receivers.wifi import dsss as wifi_dsss


def _hexbits(h):
    return np.array([int(b) for b in bin(int(h, 16))[2:].zfill(len(h) * 4)], np.uint8)


def test_adsb_crc_known_frame():
    # public example DF17 frame (ICAO 4840D6, callsign KLM1023)
    bits = _hexbits("8D4840D6202CC371C32CE0576098")
    assert adsb._residual(bits) == 0
    bits[40] ^= 1
    assert adsb._residual(bits) != 0


def _conv_encode(u):
    reg = [0, 0, 0, 0]
    out = []
    for b in u:
        d = [b] + reg
        out += [d[0] ^ d[3] ^ d[4], d[0] ^ d[1] ^ d[3] ^ d[4]]
        reg = [b] + reg[:3]
    return np.array(out, np.uint8)


def _crc10(bits25):
    poly = 0b10101110101
    reg = list(bits25) + [0] * 10
    for i in range(25):
        if reg[i]:
            for j in range(11):
                reg[i + j] ^= (poly >> (10 - j)) & 1
    return np.array(reg[25:], np.uint8) ^ 1


def test_gsm_sch_coding_roundtrip():
    rng = np.random.default_rng(1)
    d = rng.integers(0, 2, 25).astype(np.uint8)
    u = np.r_[d, _crc10(d), np.zeros(4, np.uint8)]
    c = _conv_encode(u)
    soft = 1.0 - 2.0 * c
    soft[[3, 20, 41, 60]] *= -1                       # 4 channel errors
    out = gsm.viterbi(soft)
    assert np.array_equal(out[:25], d)
    assert gsm.crc10_ok(out[:25], out[25:35])


def test_dsss_descrambler_inverts_scrambler():
    rng = np.random.default_rng(2)
    data = rng.integers(0, 2, 400).astype(np.uint8)
    reg, tx = [1, 1, 0, 1, 1, 0, 0], []
    for b in data:
        s = b ^ reg[3] ^ reg[6]
        tx.append(s)
        reg = [s] + reg[:6]
    rx = wifi_dsss._descramble(np.array(tx, np.uint8))
    assert np.array_equal(rx[7:], data[7:])            # self-synchronising after 7 bits


def test_clipping_counts_and_normalised():
    from spectrum_intel.detect.artifacts import clipping
    rng = np.random.default_rng(1)
    x = (rng.normal(0, 200, 10000) + 1j * rng.normal(0, 200, 10000)).astype(np.complex64)
    assert clipping(x)["clipped_fraction"] == 0.0                       # raw counts, no clipping
    x[:10] = 2047 + 0j
    assert clipping(x)["clipped_fraction"] == 10 / len(x)
    y = x / 2048                                                        # the same block normalised to [-1, 1)
    y[10:] = (rng.normal(0, 0.1, len(y) - 10) + 1j * rng.normal(0, 0.1, len(y) - 10)).astype(np.complex64)
    assert clipping(y)["clipped_fraction"] == 10 / len(y)               # (was 1.0: threshold 0.995 x (1 - 1))


def test_gsm_fcch_with_unknown_cfo():
    rng = np.random.default_rng(2)
    n = int(0.03 * gsm.OFS)
    y = (rng.normal(0, 1, n) + 1j * rng.normal(0, 1, n)) * 0.3
    s, cfo = 9000, 9150.0                                               # +9.7 ppm at 943 MHz
    k = np.arange(148 * gsm.SPS)
    y[s:s + len(k)] += np.exp(2j * np.pi * (gsm.F_FCCH + cfo) * k / gsm.OFS)
    got = gsm.find_fcch_cfo(y, 30e-6 * 943e6 + 2e3)
    assert got is not None
    s0, c, pur = got
    assert abs(s0 - s) <= 6 * gsm.SPS and abs(c - cfo) < 30 and pur > 0.8   # SCH search spans +-20 symbols
    assert gsm.find_fcch(y) is None                                     # without the CFO search the tone is missed
    assert gsm.find_fcch_cfo(y[:s], 30e-6 * 943e6 + 2e3) is None        # noise only


def _radar(fs=15.36e6, n=1 << 19, off=1.0e6, pris=(1120e-6, 1150e-6, 1170e-6), sweep=0.9e6, tau=69e-6, seed=3):
    """Pulse groups on a staggered PRI: 1 us short pulse, then a tau LFM chirp (S-band ATC style), in noise."""
    rng = np.random.default_rng(seed)
    x = ((rng.normal(0, 1, n) + 1j * rng.normal(0, 1, n)) * 0.05).astype(np.complex64)
    t0, k = 200e-6, 0
    while t0 + 80e-6 < n / fs:
        a = int(t0 * fs)
        x[a:a + int(1e-6 * fs)] += np.exp(2j * np.pi * off * np.arange(int(1e-6 * fs)) / fs)
        b = a + int(5e-6 * fs)
        tt = np.arange(int(tau * fs)) / fs
        x[b:b + len(tt)] += np.exp(2j * np.pi * ((off - sweep / 2) * tt + sweep / (2 * tau) * tt ** 2))
        t0 += pris[k % len(pris)]; k += 1
    return x, fs


def test_pulse_train_staggered_lfm_radar():
    from spectrum_intel.classify import rules
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.features.measure import pulse_train
    from spectrum_intel.model import Signal
    x, fs = _radar()
    fc = 2887.57e6
    s = Signal(bursts=[], f0=fc + 0.2e6, f1=fc + 1.8e6, snr_db=20.0)        # 2888.57 MHz +- 0.8 MHz
    f = pulse_train(s, Channelizer(x, fs, fc))
    assert f["pulse_n"] >= 25 and 60e-6 <= f["pulse_width_s_median"] <= 80e-6
    assert abs(f["pri_s"] - 1150e-6) < 30e-6 and f["pri_spread"] > 0.02          # staggered
    assert abs(f["lfm_sweep_hz"] - 0.9e6 * 0.8) < 0.15e6 and f["lfm_r2"] > 0.95   # central 80 % of the chirp
    f.update(continuous=False, n_bursts=3, burst_dur_s_median=320e-6)            # spectrogram boxes: too long
    name, sc, why = rules.r_pulsed(s, f)      # long LFM pulse: solid-state ATC (research, signal gaps S10)
    assert name == "Pulsed radar (S-band ATC, solid-state)-like" and sc == 0.8 and any("LFM" in w for w in why)


def test_narrow_line_has_no_pulse_train():
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.features.measure import pulse_train
    from spectrum_intel.model import Signal
    x, fs = _radar()
    s = Signal(bursts=[], f0=2887.57e6 + 1.0e6, f1=2887.57e6 + 1.0188e6)      # 18.8 kHz: no us pulse fits
    assert pulse_train(s, Channelizer(x, fs, 2887.57e6)) == {}


def test_bluetooth_label_by_channel():
    from spectrum_intel.core.emission import Emission
    from spectrum_intel.receivers.rc_link.plugin import RcLinkDecoder
    em = lambda c: [Emission("bluetooth", "bluetooth_like", 0, 300, c - 0.5, c + 0.5, c)]
    d = RcLinkDecoder()
    assert d.pretty("bluetooth_like", em(2426.02)) == "Bluetooth LE advertising-like"     # channel 38
    assert d.pretty("bluetooth_like", em(2466.97)) == "Bluetooth (BLE/BR)-like"           # data channel / BR hop


def test_pulsed_rule_gates():
    """Pulse-train gates (live 1090 MHz captures): noise peaks too short for the band, single pulses in the DME band
    and anything on the SSR frequencies are not DME; a 15 us pulse pair at 1000 MHz is."""
    from spectrum_intel.classify import rules
    from spectrum_intel.model import Signal
    f = dict(continuous=False, n_bursts=8, burst_dur_s_median=2e-3, pulse_n=700, pulse_width_s_median=2e-6)
    assert rules.r_pulsed(Signal(bursts=[], f0=999.8e6, f1=1000.2e6), f) is None          # 0.4 MHz x 2 us < 1
    f.update(pulse_width_s_median=6e-6)
    assert rules.r_pulsed(Signal(bursts=[], f0=999.5e6, f1=1000.5e6), f) is None          # one pulse, not a pair
    f.update(pulse_width_s_median=15.5e-6)
    assert rules.r_pulsed(Signal(bursts=[], f0=1089.6e6, f1=1090.4e6), f) is None         # SSR reply frequency
    assert rules.r_pulsed(Signal(bursts=[], f0=999.6e6, f1=1000.4e6), f)[0] == "DME/TACAN pulse pairs-like"


def test_gsm_burst_is_not_an_rc_link():
    """Live LO 940 (LIVE_GSM_958MHz_fsk): GSM-900 traffic bursts (554 us GMSK on ARFCN 110 / 115) were labelled
    sik-like. One to eight slots on the 200 kHz raster in a GSM band are GSM; off-raster or outside GSM, not."""
    from spectrum_intel.receivers.rc_link.detector import _gsm_burst
    r = dict(t0_us=10058.4, t1_us=10612.4)
    assert _gsm_burst(r, 958.0079e6)                         # one burst, +7.9 kHz (receiver ppm)
    assert _gsm_burst(dict(t0_us=0.0, t1_us=3 * 576.9), 945.2e6)   # three adjacent slots
    assert not _gsm_burst(r, 958.1e6)                        # between two ARFCNs
    assert not _gsm_burst(r, 866.0e6)                        # 865-867 MHz SRD: an RC / telemetry link may sit there
    assert not _gsm_burst(dict(t0_us=0.0, t1_us=2000.0), 958.0e6)  # not a whole number of slots


def test_pulse_train_pieces_merge():
    """Live 3036-3047 MHz emitter (LIVE_PULSED_3040MHz): three narrower detections share one 4 ms pulse train ->
    one signal; a sliver on its pulses joins, a burst off the pulses and another PRI stay apart."""
    from spectrum_intel.detect.merge import merge_pulse_trains
    from spectrum_intel.model import Box, Signal
    t = [1.548e-3 + 4e-3 * k for k in range(9)]

    def sig(f0, f1, pulses=None, pri=None, bursts=None):
        s = Signal(bursts=bursts or [Box(p, p + 10e-6, f0, f1, 5, 3, "fine") for p in (pulses or t)], f0=f0, f1=f1,
                   snr_db=3.0, detection={"description": "x"})
        if pri:
            s.features.update(pri_s=pri, _pulse_t_s=list(pulses or t))
        return s
    a, b = sig(3036.8e6, 3039.5e6, pri=4e-3), sig(3040.1e6, 3041.8e6, [v + 2e-6 for v in t], pri=4e-3)
    sliver = sig(3036.7e6, 3037.1e6, bursts=[Box(t[3] - 5e-6, t[3] + 12e-6, 3036.7e6, 3037.1e6, 4, 2, "fine")])
    off = sig(3038.4e6, 3038.7e6, bursts=[Box(0.0, 83e-6, 3038.4e6, 3038.7e6, 4, 2, "fine")])
    other = sig(3029e6, 3030e6, [v + 1e-3 for v in t], pri=4e-3)
    out, merged = merge_pulse_trains([a, off, b, sliver, other])
    assert len(out) == 3 and len(merged) == 1
    m = merged[0]
    assert (m.f0, m.f1) == (3036.7e6, 3041.8e6) and len(m.bursts) == 19
    assert "3 energy detections with one pulse train" in m.detection["description"]
    assert off in out and other in out


def test_short_pulses_are_not_ofdm():
    """An autocorrelation lag longer than a quarter of the pulses cannot be a CP-OFDM symbol."""
    from spectrum_intel.features.measure import infer_modulation
    from spectrum_intel.model import Signal
    s = Signal(bursts=[], f0=3036.8e6, f1=3041.8e6, snr_db=8.0)
    f = dict(ofdm_scs_hz=180410.0, acf_peak=0.08, acf_peak_over_median=17.1, acf_peak_lag_s=5.54e-6)
    assert infer_modulation(s, dict(f, pulse_n=9, pulse_width_s_median=11e-6), 1e3)[0] != "OFDM"
    assert infer_modulation(s, dict(f, pulse_n=9, pulse_width_s_median=200e-6), 1e3)[0] == "OFDM"
    assert infer_modulation(s, f, 1e3)[0] == "OFDM"


def test_wifi_offgrid_link_found():
    """Live LO 4877 (LIVE_MIMOSA_4880MHz): an 802.11 link at 4880 MHz, off the channel grid and outside the Wi-Fi
    bands. The CP trigger plans blind centres around it and the detector decodes its PPDUs; on a standard channel
    the decoded PPDUs are blanked and nothing is planned; noise plans nothing."""
    from spectrum_intel.receivers.wifi import ofdm_rx as W, ofdm_tx as TX
    ack = bytes([0xD4, 0, 0, 0]) + bytes.fromhex("20b5c646e83e")
    burst = TX.ppdu(ack, fmt="NON_HT", rate_mbps=6)
    fs, n = 61.44e6, 1 << 19

    def scene(f_hz, fc):
        return TX.scene([dict(samples=burst, f_hz=f_hz, t_s=50e-6 + 450e-6 * k, snr_db=15) for k in range(18)],
                        n, fs, fc, seed=4)
    x = scene(4880e6, 4877e6)
    plan = W.offgrid_plan(x, fs, 4877e6, [])
    assert any(abs(c["freq_hz"] - 4880e6) < 1 for c in plan) and not any(c["standard"] for c in plan)
    pk = W.merge_packets([], W.WifiDetector(fs, 4877e6, channels=plan).run(x))
    assert len(pk) >= 15 and all(abs(p.center_mhz - 4880.0) < 0.1 for p in pk)
    xs = scene(5180e6, 5190e6)                              # standard channel 36: the normal scan owns it
    known = W.WifiDetector(fs, 5190e6).run(xs)
    assert len(known) >= 15 and W.offgrid_plan(xs, fs, 5190e6, known) == []
    rng = np.random.default_rng(1)
    noise = ((rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)).astype(np.complex64)
    assert W.offgrid_plan(noise, fs, 4877e6, []) == []


def test_wifi_offgrid_weak_low_duty():
    """PLAN G13 (LIVE_MIMOSA_4880MHz_weak): a weak, low-duty off-grid link is diluted in the whole-capture CP
    correlation; the windowed fallback (250 us windows, gated by a whole-capture ratio >= 2) still plans it, and noise
    plans nothing."""
    from spectrum_intel.receivers.wifi import ofdm_rx as W, ofdm_tx as TX
    fs, n = 61.44e6, 1 << 21
    burst = TX.ppdu(bytes([0x08, 0x01]) + bytes(1000), fmt="NON_HT", rate_mbps=6)
    x = TX.scene([dict(samples=burst, f_hz=4880e6, t_s=1e-3 + 10e-3 * k, snr_db=0) for k in range(3)], n, fs, 4877e6,
                 seed=4)
    win = W._offgrid_windows
    try:
        W._offgrid_windows = lambda *a: False                 # the whole-capture test alone misses it
        assert not any(abs(c["freq_hz"] - 4880e6) < 1 for c in W.offgrid_plan(x, fs, 4877e6, []))
    finally:
        W._offgrid_windows = win
    assert any(abs(c["freq_hz"] - 4880e6) < 1 for c in W.offgrid_plan(x, fs, 4877e6, []))
    rng = np.random.default_rng(2)
    noise = ((rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)).astype(np.complex64)
    assert W.offgrid_plan(noise, fs, 4877e6, []) == []


def test_coincident_pieces_join_their_receiver():
    """Live LO 5755 (LIVE_WIFI_5755MHz): a wideband 802.11 PPDU split into pieces at its fading notches. >= 3 pieces in
    the Wi-Fi signal's span with the same bursts and >= 8 MHz together are absorbed; a lone narrow hop and pieces at
    other times are not."""
    from spectrum_intel.analyze import _absorb_coincident
    from spectrum_intel.model import Box, Signal

    def sig(f0, f1, t=((0.010, 0.0106),)):
        return Signal(bursts=[Box(a, b, f0, f1, 20, 15, "fine") for a, b in t], f0=f0, f1=f1)
    wifi = sig(5736.7e6, 5753.3e6, ((0.0101, 0.0102), (0.020, 0.0201)))
    pieces = [sig(5737e6, 5740e6), sig(5741e6, 5743e6, ((0.010005, 0.010598),)), sig(5746e6, 5747e6), sig(5749e6, 5752e6)]
    hop = sig(5744e6, 5745e6, ((0.0150, 0.0152),))
    late = [sig(5737e6, 5739e6, ((0.0300, 0.0306),)), sig(5740e6, 5741e6, ((0.0300, 0.0306),))]
    taken = {}
    _absorb_coincident([wifi], pieces + [hop] + late, taken)
    assert [id(e) for e in taken[0]] == [id(e) for e in pieces]


def test_no_radar_label_on_wifi_fragments():
    """5495.65 MHz (U-NII-2e, shared with Wi-Fi): bursty fragments without a PRI are not radar; with a PRI they are."""
    from spectrum_intel.classify import rules
    from spectrum_intel.model import Signal
    s = Signal(bursts=[], f0=5495.5e6, f1=5495.8e6)
    f = dict(continuous=False, n_bursts=4, burst_dur_s_median=20e-6)
    assert rules.r_pulsed(s, f) is None
    assert rules.r_pulsed(s, dict(f, burst_period_s=1e-3))[0] == "Pulsed radar (C-band)-like"   # S10 class


def _esb_frames(n, plen=8, seed=3):
    """Standard nRF ESB frames as bits from the address on: 5-byte address (MSB '11', so the sync word starts at the
    address after a 1010.. preamble), 9-bit PCF, payload with a counter, CRC-16/CCITT init 0xFFFF over all of it."""
    rng = np.random.default_rng(seed)
    addr = np.r_[1, 1, rng.integers(0, 2, 38)].astype(np.uint8)
    out = []
    for k in range(n):
        pcf = np.array([(plen >> (5 - i)) & 1 for i in range(6)] + [(k >> 1) & 1, k & 1, 0], np.uint8)
        pay = np.r_[[(k >> (7 - i)) & 1 for i in range(8)], rng.integers(0, 2, 8 * plen - 8)].astype(np.uint8)
        m = np.r_[addr, pcf, pay]
        c = 0xFFFF
        for b in m:
            fb = ((c >> 15) & 1) ^ int(b)
            c = ((c << 1) & 0xFFFF) ^ (0x1021 if fb else 0)
        out.append(np.r_[m, [(c >> (15 - i)) & 1 for i in range(16)]].astype(np.uint8))
    return out


def test_esb_crc_found_without_knowing_the_frame():
    """esb.verify_links: the CRC end and init of standard ESB frames (init 0xFFFF at the address) come out of the
    packets alone; a corrupted packet fails; frames whose CRC fields are random never pass."""
    from spectrum_intel.receivers.rc_link import esb
    fr = _esb_frames(12)
    fr[5] = fr[5].copy(); fr[5][60] ^= 1
    v = esb.verify_links(fr, [1000.0 * k for k in range(12)])
    assert [d["crc_ok"] for d in v] == [k != 5 for k in range(12)]
    assert v[0]["crc_span_bits"] == 40 + 9 + 64 and v[0]["crc_init_hex"] == "FFFF"
    assert abs(v[0]["node_frame_period_us"] - 1000.0) < 1e-6
    rng = np.random.default_rng(7)
    bad = [np.r_[f[:-16], rng.integers(0, 2, 16)].astype(np.uint8) for f in fr]
    assert not any(d["crc_ok"] for d in esb.verify_links(bad, [1000.0 * k for k in range(12)]))


def test_esb_link_detected_and_crc_checked():
    """A 1 kHz 2 Mb/s GFSK ESB link (h 0.5, 16-bit preamble) on one channel, 2446 MHz, for a whole block at 15 dB:
    the packets are ESB-like bursts on the 1 MHz raster and pass the joint CRC check. (Before, fsk_refine took a
    ~70 us packet in fsk_scan's ~870 us box for noise, and the 16-bit preamble measured 14 bits.) 12 dB: 7 of 16
    pass, 9 dB: none."""
    import scipy.fft as sfft
    from spectrum_intel.receivers.rc_link import detector as D
    from spectrum_intel.receivers.rc_link.synth import gfsk_burst
    fs, fc, n = 61.44e6, 2440e6, 1 << 21
    rng = np.random.default_rng(11)
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)
    t = np.arange(n) / fs
    amp = np.sqrt(10 ** 1.5 * 2e6 / fs)                    # 15 dB in a 2 MHz band
    fr = _esb_frames(16)
    for k, f in enumerate(fr):
        b, _ = gfsk_burst(2e6, 500e3, bits=np.r_[np.arange(16) % 2 == 0, f, 0, 1, 0, 1].astype(np.uint8), fs=fs, rng=rng)
        i = int((0.5e-3 + k * 1e-3) * fs)
        x[i:i + len(b)] += amp * b * np.exp(2j * np.pi * 6e6 * t[i:i + len(b)]).astype(np.complex64)
    ems = [e for e in D.detect_fsk(x, sfft.fft(x), fs, fc) if e.family == "esb_link"]
    assert len(ems) >= 15 and all(round(e.center_mhz) == 2446 for e in ems)
    assert sum(e.decoded["crc_ok"] for e in ems) >= 14
    d = next(e.decoded for e in ems if e.decoded["crc_ok"])
    assert d["crc_span_bits"] == 113 and abs(d["node_frame_period_us"] - 1000.0) < 1.0
