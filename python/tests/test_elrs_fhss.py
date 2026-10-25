"""ExpressLRS FHSS hop sequence (research, PLAN R2d; research/drone_transmitters/07_elrs_details.md section 4):
sequence structure, UID2 / UID3 recovery from observed hops, hop recovery from a SYNC packet."""
import random

from spectrum_intel.receivers.rc_link import elrs as E


def test_sequence_structure():
    for count in (80, 40, 13, 4):
        seq = E.fhss_sequence(E.fhss_seed(1, 2, 3, 4), count)
        n = (256 // count) * count
        assert len(seq) == n
        for b in range(0, n, count):                       # every block: sync first, a permutation of all channels
            assert seq[b] == count // 2 and sorted(seq[b:b + count]) == list(range(count))
    a = E.fhss_sequence(E.fhss_seed(0x12, 2, 3, 4), 40)
    assert a == E.fhss_sequence(E.fhss_seed(0x92, 2, 3, 4), 40)     # seed bit 31 is dropped by the LCG
    assert a != E.fhss_sequence(E.fhss_seed(0x13, 2, 3, 4), 40)


def test_uid23_recovery():
    rnd = random.Random(3)
    for count in (40, 13):
        uid = [rnd.randrange(256) for _ in range(6)]
        seq = E.fhss_sequence(E.fhss_seed(*uid[2:6]), count)
        k0 = rnd.randrange(256)
        obs = [(k0 + d, seq[(k0 + d) % len(seq)]) for d in range(6)]
        c = E.uid23_candidates(uid[4], uid[5], count, obs)
        assert (uid[2], uid[3]) in c and len(c) <= 32
        assert (uid[2], uid[3]) in E.uid23_candidates(uid[4], uid[5], count, obs, uid3=uid[3])


def test_recover_hops_from_sync():
    rnd = random.Random(5)
    uid = [rnd.randrange(256) for _ in range(6)]
    count, hop, T = 40, 2, 5000.0
    seq = E.fhss_sequence(E.fhss_seed(*uid[2:6]), count)
    n0, i0 = 17, 33
    idx = lambda dn: (i0 + sum(1 for n in range(n0, n0 + dn) if (n + 1) % hop == 0)) % 240
    pk = [(1000 + T * dn, seq[idx(dn)]) for dn in range(1, 13)]
    r = E.recover_hops(dict(sync=dict(fhss_index=i0, nonce=n0, uid4=uid[4], uid5=uid[5]), ota_version=4),
                       1000.0, seq[i0], T, hop, count, pk)
    assert f"{uid[2]:02X}{uid[3]:02X}" in r["uid2_uid3"] or r["uid23_candidates"] > 8
    assert r["observations"] == 13


def test_hop_recovery_in_capture():
    """A SYNC packet + the link's following packets on the channels its UID's sequence gives -> hop_recovery on the
    SYNC emission holds the true UID2 / UID3 among its candidates."""
    import numpy as np
    from spectrum_intel.receivers.rc_link import detector as D, synth as S
    from spectrum_intel.receivers.rc_link.plugin import PLUGINS
    uid = [1, 2, 0x5A, 0x3C, 0x77, 0x21]
    band, sf, bw, cr_den, pre, pl, hop, ival = D.ELRS_LORA["900_100hz"]
    fs, fc, count = 30.72e6, 915.2e6, 40
    seq = E.fhss_sequence(E.fhss_seed(*uid[2:6]), count)
    rng = np.random.default_rng(4)
    x = S.noise_capture(1 << 21, rng)
    n0, i0 = 11, 120
    i0 -= i0 % count                                       # a block start: the sync channel
    bursts = []
    for dn in range(7):
        idx = (i0 + sum(1 for n in range(n0, n0 + dn) if (n + 1) % hop == 0)) % len(seq)
        pkt = E.build_sync(i0, n0, 2, 3, uid) if dn == 0 else bytes(rng.integers(0, 256, pl).tolist())
        iq, _ = S.lora_packet(sf, bw, n_pre=pre, payload=S.lora_implicit_symbols(pkt, sf, cr_den - 4), sync_word=0x12,
                              cfo_hz=3e3, fs=fs, rng=rng)
        f = 903.5e6 + (926.9e6 - 903.5e6) * seq[idx] / (count - 1)
        bursts.append(dict(t0_us=1000 + ival * dn, f_hz=f, iq=iq, bw_hz=bw))
    y, _ = S.embed(x, fs, fc, bursts, 12)
    em = D.detect(y, fs, fc)
    syn = [e for e in em if (e.decoded or {}).get("packet_type") == "SYNC"]
    assert syn and syn[0].decoded["crc_ok"]
    r = syn[0].decoded["hop_recovery"]
    assert r["domain"] == "FCC915" and r["observations"] >= 5
    assert r["uid23_candidates"] > 8 or "5A3C" in r["uid2_uid3"]
    f = PLUGINS[0].summarise([e for e in em if "elrs" in e.label], None)["fields"]["elrs_sync"]
    assert f["hop_recovery"]["uid23_candidates"] == r["uid23_candidates"]
