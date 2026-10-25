"""DJI DroneID frame variants (research, PLAN D4; research/drone_transmitters/09_ble_coded_droneid_legacy.md):
9-symbol (OcuSync 2/3) vs 8-symbol legacy (Mavic 2 / older) on synthetic bursts through dji detect."""
import numpy as np
import pytest

from spectrum_intel.receivers.dji import detector as D, synth as DS
from spectrum_intel.receivers.rc_link import synth as RS

FS, FC = 61.44e6, 2440e6


def run(ns, snr_db, f_hz, seed):
    rng = np.random.default_rng(seed)
    b = DS.upsample(DS.droneid_burst(ns, rng), 4)
    x = RS.noise_capture(1 << 21, rng)
    y, _ = RS.embed(x, FS, FC, [dict(t0_us=5000, f_hz=f_hz, iq=b, bw_hz=9e6)], snr_db)
    return [e for e in D.detect(y, FS, FC, ctx={"decode_droneid": False}) if e.family == "dji_droneid"]


@pytest.mark.parametrize("ns", [9, 8])
@pytest.mark.parametrize("snr_db,f_hz,seed", [(10, 2444.5e6, 1), (3, 2431.0e6, 2), (10, 2440.0e6, 3)])
def test_droneid_frame_variant(ns, snr_db, f_hz, seed):
    em = run(ns, snr_db, f_hz, seed)
    assert len(em) == 1, [(e.start_us, e.center_mhz) for e in em]
    e = em[0]
    assert e.features["n_symbols"] == ns
    assert abs(e.start_us - 5000) < 2 and abs(e.center_mhz - f_hz / 1e6) < 0.02
    assert abs(e.features["burst_us"] - (643.2 if ns == 9 else 571.4)) < 2
