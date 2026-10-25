"""Cross-check elrs.fhss_sequence against a line-by-line transliteration of the GPL ExpressLRS FHSS.cpp / random.cpp
(oracle only: this bench script is not part of the library)."""
import random
import sys
sys.path.insert(0, "python")
from spectrum_intel.receivers.rc_link import elrs as E

def oracle(seed, freqCount):
    st = [seed]
    def rng():
        st[0] = (214013 * st[0] + 2531011) % 2147483648
        return st[0] >> 16
    seqCount = (256 // freqCount) * freqCount
    sync = freqCount // 2
    seq = [0] * seqCount
    for i in range(seqCount):
        seq[i] = sync if i % freqCount == 0 else (0 if i % freqCount == sync else i % freqCount)
    for i in range(seqCount):
        if i % freqCount != 0:
            offset = ((i // freqCount) * freqCount) & 0xFF
            r = (rng() % (freqCount - 1) + 1) & 0xFF
            seq[i], seq[(offset + r) & 0xFF] = seq[(offset + r) & 0xFF], seq[i]
    return seq

rnd = random.Random(1)
bad = 0
for count in (80, 40, 20, 13, 8, 4, 3):
    for _ in range(50):
        seed = rnd.getrandbits(32)
        bad += E.fhss_sequence(seed, count) != oracle(seed, count)
print("mismatches", bad)
