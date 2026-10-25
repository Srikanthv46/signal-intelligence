"""TBS Crossfire 150 Hz FSK (research, PLAN R9; research/drone_transmitters/08_crossfire_dsm.md): synthetic uplink
packets on the 868 grid through rc_link detect -> sync word, packet bytes, grid channel."""
import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import crossfire as CF, detector as D, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 868e6


def test_grid():
    assert CF.grid_channel(860.165e6 + 7 * 260e3) == {"band": "868", "channel": 7, "link": "uplink"}
    assert CF.grid_channel(860.165e6 + 57 * 260e3) == {"band": "868", "channel": 7, "link": "downlink"}
    assert CF.grid_channel(861.0e6) is None


@pytest.mark.parametrize("snr_db", [20, 12])
def test_crossfire_uplink(snr_db):
    rng = np.random.default_rng(5)
    x = S.noise_capture(1 << 21, rng)
    data = bytes(rng.integers(0, 256, 23).tolist())
    bursts = []
    for k, n in enumerate((12, 30, 44)):                      # three hops, 6.666 ms apart
        iq, _ = S.crossfire_burst(data, cfo_hz=3e3, rng=rng)
        bursts.append(dict(t0_us=1000 + 6666 * k, f_hz=860.165e6 + n * 260e3, iq=iq, bw_hz=170e3))
    y, _ = S.embed(x, FS, FC, bursts, snr_db)
    em = [e for e in D.detect(y, FS, FC) if e.label == "crossfire_sync"]
    assert len(em) == 3
    for e in em:
        d = e.decoded
        assert d["sync"] == "0x3E6B3E" and d["grid"]["link"] == "uplink"
        errs = bin(int(d["dewhitened_hex"], 16) ^ int(data.hex(), 16)).count("1")
        assert errs <= (0 if snr_db >= 20 else 3), errs              # hard slicer, no CRC to correct with
    f = PLUGINS[0].summarise(em, None)["fields"]["crossfire"]
    assert f["packets"] == 3 and f["band"] == "868" and [c for _, c in f["channels"]] == [12, 30, 44]
    assert PLUGINS[0].pretty("crossfire_sync", em).startswith("TBS Crossfire-like")
