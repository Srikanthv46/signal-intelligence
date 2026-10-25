"""Blind FSK frame structure (research, PLAN F1): preamble length, 32-bit sync candidate, payload length and known
sync-word naming on synthetic packets, direct and through rc_link detect + the plugin's fsk_frame field."""
import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import detector as D, fsk_frame as FF, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 2440e6


def _bits(sync, nsync, n_pay, rng, pre=32):
    return np.r_[np.arange(pre) % 2 == 0, [(sync >> (nsync - 1 - i)) & 1 for i in range(nsync)],
                 rng.integers(0, 2, n_pay), np.zeros(8)].astype(np.uint8)


@pytest.mark.parametrize("sync,n,name", [(0x5475C52A, 32, "FlySky AFHDS 2A-like (A7105 ID code)"),
                                         (0x9A3C5E71, 32, None)])
def test_frame_direct(sync, n, name):
    rng = np.random.default_rng(2)
    bits = _bits(sync, n, 200, rng)
    iq, _ = S.gfsk_burst(250e3, 125e3, bits=bits, fs=2e6, rng=rng)
    for y in (iq, np.conj(iq)):
        fr = FF.frame(y, 2e6, 250e3)
        assert fr["preamble_bits"] >= 28 and abs(fr["payload_bits"] - 208) <= 34
        assert (fr.get("known_sync") or {}).get("name") == name
    if name is None:   # the alternating run absorbs the sync's leading "10": the candidate starts 2 bits in
        assert FF.frame(iq, 2e6, 250e3)["sync32_hex"][:6] == f"{(0x9A3C5E71 << 2) & 0xFFFFFFFF:08X}"[:6]


def test_frame_detect_and_field():
    rng = np.random.default_rng(4)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k, f in enumerate((2425.3e6, 2433.1e6, 2451.9e6)):
        iq, _ = S.gfsk_burst(250e3, 125e3, bits=_bits(0x9A3C5E71, 32, 240, rng), rng=rng)
        bursts.append(dict(t0_us=2000 + 5000 * k, f_hz=f, iq=iq, bw_hz=600e3))
    iq, _ = S.gfsk_burst(250e3, 125e3, bits=_bits(0xD391D391, 32, 240, rng), rng=rng)     # CC2500 sync, not FrSky rate
    bursts.append(dict(t0_us=20000, f_hz=2445.0e6, iq=iq, bw_hz=600e3))
    y, _ = S.embed(x, FS, FC, bursts, 18)
    em = [e for e in D.detect(y, FS, FC) if e.features.get("frame")]
    assert len(em) == 4
    f = PLUGINS[0].summarise(em, None)["fields"]["fsk_frame"]
    assert f["packets"] == 4 and f["sync24_candidates"][0] == ["68F179", 3]
    assert f["known_sync"] == ["FrSky D16-like (CC2500)"]
