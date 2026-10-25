"""ExpressLRS identity bit (research, PLAN R2a; research/drone_transmitters/07_elrs_details.md): on 2.4 GHz the SX1280
inverts IQ when UID[5] bit 0 = 1 -> plugin field elrs_identity."""
import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import detector as D, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 2440e6


@pytest.mark.parametrize("inv", [True, False])
def test_uid5_bit0_from_iq_polarity(inv):
    rng = np.random.default_rng(11 if inv else 12)
    x = S.noise_capture(1 << 21, rng)
    b = S.elrs_link("2g4_500hz", rng=rng, iq_invert=inv, f_window=(FC - 26e6, FC + 26e6), dur_s=len(x) / FS - 1e-3)
    y, _ = S.embed(x, FS, FC, b, 15)
    em = [e for e in D.detect(y, FS, FC) if "elrs" in e.label]
    assert em
    f = PLUGINS[0].summarise(em, None)["fields"]
    assert f["elrs_identity"]["uid5_bit0"] == int(inv) and f["elrs_identity"]["packets"] == len(em)
