#!/usr/bin/env python3
"""Reference for the ExpressLRS FHSS functions (rc_link/elrs; research, PLAN R2d): sequences, UID2 / UID3 candidate
searches and hop recovery from a SYNC packet.

    .venv/bin/python -I c++/tests/vectors/elrs_fhss/gen_elrs_fhss.py
"""
import json
import random
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import elrs as E   # noqa: E402

rnd = random.Random(91)
seqs = [[s, c, E.fhss_sequence(s, c)] for s, c in ((0, 80), (0xFFFFFFFF, 40), (0x12345678, 13), (rnd.getrandbits(32), 4),
                                                    (rnd.getrandbits(32), 3), (rnd.getrandbits(32), 20))]
searches, recs = [], []
for count, n_obs, u3k in ((40, 6, False), (13, 6, False), (8, 4, True), (40, 3, True)):
    uid = [rnd.randrange(256) for _ in range(6)]
    seq = E.fhss_sequence(E.fhss_seed(*uid[2:6]), count)
    k0 = rnd.randrange(256)
    obs = [[k0 + d, seq[(k0 + d) % len(seq)]] for d in range(n_obs)]
    args = [uid[4], uid[5], count, obs, 4, uid[3] if u3k else None]
    searches.append([args, E.uid23_candidates(*args)])
for count, hop, T, n in ((40, 2, 5000.0, 12), (13, 4, 4000.0, 8), (4, 2, 20000.0, 3)):
    uid = [rnd.randrange(256) for _ in range(6)]
    seq = E.fhss_sequence(E.fhss_seed(*uid[2:6]), count)
    n0, i0 = rnd.randrange(256), rnd.randrange(len(seq))
    ln = len(seq)
    idx = lambda dn: (i0 + sum(1 for m in range(n0, n0 + dn) if (m + 1) % hop == 0)) % ln
    pk = [[1000 + T * dn + 30, seq[idx(dn)]] for dn in range(1, n + 1)]
    sync = dict(sync=dict(fhss_index=i0, nonce=n0, uid4=uid[4], uid5=uid[5]), ota_version=4)
    recs.append([[sync, 1000.0, seq[i0], T, hop, count, pk], E.recover_hops(sync, 1000.0, seq[i0], T, hop, count, pk)])
(OUT / "cases.json").write_text(json.dumps(dict(seqs=seqs, searches=searches, recover=recs)))
print([len(s[1]) for s in searches], [r[1] and r[1]["uid23_candidates"] for r in recs])

# whole rc_link pipeline (PLAN R2d): IN866 (4 channels), 900 100 Hz mode, a SYNC packet + the following packets on the
# channels the UID's sequence gives; 7.68 MSps, 2^18 samples
import numpy as np                                                      # noqa: E402
from spectrum_intel.receivers.rc_link import detector as D, synth as S   # noqa: E402
uid = [1, 2, 0x5A, 0x3C, 0x77, 0x21]
band, sf, bw, cr_den, pre, pl, hop, ival = D.ELRS_LORA["900_100hz"]
fs, fc, count = 7.68e6, 866.16e6, 4
seq = E.fhss_sequence(E.fhss_seed(*uid[2:6]), count)
rng = np.random.default_rng(93)
x = S.noise_capture(1 << 18, rng)
n0, i0 = 11, 96
bursts = []
for dn in range(3):
    idx = (i0 + sum(1 for n in range(n0, n0 + dn) if (n + 1) % hop == 0)) % len(seq)
    pkt = E.build_sync(i0, n0, 2, 3, uid) if dn == 0 else bytes(rng.integers(0, 256, pl).tolist())
    iq, _ = S.lora_packet(sf, bw, n_pre=pre, payload=S.lora_implicit_symbols(pkt, sf, cr_den - 4), sync_word=0x12,
                          cfo_hz=2e3, fs=fs, rng=rng)
    bursts.append(dict(t0_us=1000 + ival * dn, f_hz=865.375e6 + (866.95e6 - 865.375e6) * seq[idx] / 3, iq=iq, bw_hz=bw))
y, _ = S.embed(x, fs, fc, bursts, 12)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, fs, fc) if (e.decoded or {}).get("packet_type") == "SYNC"]
ref = json.loads((OUT / "cases.json").read_text())
ref["pipeline"] = dict(fs=fs, fc=fc, sync=[e.decoded for e in em])
(OUT / "cases.json").write_text(json.dumps(ref))
print("pipeline", [e.decoded.get("hop_recovery") for e in em])
