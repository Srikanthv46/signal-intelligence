"""Reference vectors for c++/tests/unit/test_rx_nr.cpp, written by the Python receiver
(spectrum_intel.receivers.nr.ssb, spectrum_intel.core.seq3gpp) and by py3gpp used as a
black-box oracle only (AGPL-3: nothing of it is copied into the C++ library).

    .venv/bin/python c++/tests/vectors/nr/gen_nr.py      -> nr.json next to this script
"""
import json
from pathlib import Path

import numpy as np
from py3gpp import helper
from py3gpp.nrBCH import nrBCH
from py3gpp.nrPBCHDMRS import nrPBCHDMRS
from py3gpp.nrPRBS import nrPRBS
from py3gpp.nrPolarDecode import nrPolarDecode
from py3gpp.nrPSS import nrPSS
from py3gpp.nrRateRecoverPolar import nrRateRecoverPolar
from py3gpp.nrSSS import nrSSS

from spectrum_intel.core import seq3gpp
from spectrum_intel.receivers.nr import ssb

rng = np.random.default_rng(7)
cx = lambda a: [[float(v.real), float(v.imag)] for v in np.asarray(a).ravel()]
out = {}

# Gold sequence / QPSK (core/seq3gpp) against py3gpp nrPRBS
gold = []
for ci in (0, 75, 1007, 2 ** 30 + 12345, 2 ** 31 - 1):
    g = seq3gpp.gold(ci, 3000)
    assert np.array_equal(g, np.asarray(nrPRBS(ci, 3000), int))
    gold.append({"c_init": ci, "bits": "".join(map(str, g))})
out["gold"] = gold
out["qpsk_1234"] = cx(seq3gpp.qpsk(seq3gpp.gold(1234, 288)))
out["dmrs"] = [{"pci": p, "ib": ib, "d": cx(nrPBCHDMRS(p, ib))} for p, ib in ((75, 3), (82, 0), (1007, 7))]
out["pss"] = [[int(v) for v in nrPSS(n)] for n in range(3)]
out["sss"] = [{"pci": p, "d": [int(v) for v in nrSSS(p)]} for p in (0, 82, 225, 1007)]

# polar code construction
out["reliability"] = [int(v) for v in helper.frozen_pos_table]
out["info_56_512"] = [int(v) for v in helper.generate_5g_ranking(56, 512)[1]]
out["interleaver_56"] = [int(v) for v in helper.polar_precode_interleave(56)]
llr = rng.normal(size=864)
out["rr_in"] = llr.tolist()
out["rr_out"] = nrRateRecoverPolar(llr, 56, 512, False, False).tolist()

# BCH: encoded MIBs with noise, plus pure noise, through the Python receiver's _bch_decode
bch = []
for t in range(12):
    pci = int(rng.integers(0, 1008))
    trblk = rng.integers(0, 2, 24)
    sfn, hrf = int(rng.integers(0, 1024)), int(rng.integers(0, 2))
    cw = nrBCH(trblk, sfn, hrf, 8, 0, pci)
    sigma = (0.3, 2.6, 3.5, 4.5)[t % 4]
    soft = (1 - 2 * cw) + sigma * rng.normal(size=864) if t < 10 else rng.normal(size=864)
    scr, crc, payload, lsb, h, _ = ssb._bch_decode(soft, 0, pci)
    dec = nrPolarDecode(nrRateRecoverPolar(soft, 56, 512, False, False), 56, 864, 8)
    bch.append({"pci": pci, "soft": soft.tolist(), "crc_ok": int(np.asarray(crc).ravel()[0]) == 0,
                "polar": [int(v) for v in dec], "scrblk": [int(v) for v in scr], "payload": [int(v) for v in payload],
                "sfn_lsb": [int(v) for v in lsb], "hrf": int(h)})
out["bch"] = bch

# SSB-level pieces on a random grid
G = rng.normal(size=(4, 240)) + 1j * rng.normal(size=(4, 240))
out["grid"] = cx(G)
out["sss_scores"] = {str(n2): ssb._sss_scores(G, n2).tolist() for n2 in range(3)}
llr = []
for pci in (75, 82, 3):
    ib, l = ssb._pbch_llr(G, pci)
    llr.append({"pci": pci, "ib": int(ib), "llr": l.tolist()})
out["pbch_llr"] = llr
m = rng.random(5000) ** 4
out["peaks_in"] = m.tolist()
out["peaks"] = ssb._peaks(m)
out["gscn"] = [{"lo": lo, "hi": hi, "f": ssb._gscn(lo, hi)}
               for lo, hi in ((758e6, 803e6), (2300e6, 2400e6), (3300e6, 3420e6), (2990e6, 3010e6))]
T = ssb._pss_templates(30e3, 1024)
out["tpl_30k_1024"] = cx(T)

Path(__file__).with_name("nr.json").write_text(json.dumps(out))
print("wrote nr.json")
