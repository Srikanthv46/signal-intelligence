"""The fast 3GPP sequence code (core.seq3gpp.gold / qpsk, the NR receiver's caches and its
BCH payload descrambling) gives exactly what py3gpp gives."""
import numpy as np
import pytest

py3gpp = pytest.importorskip("py3gpp")
from py3gpp.helper import _calc_gold                      # noqa: E402
from py3gpp.nrBCHDecode import nrBCHDecode                # noqa: E402
from py3gpp.nrPBCHDMRS import nrPBCHDMRS, nrPBCHDMRScinit  # noqa: E402
from py3gpp.nrPBCHPRBS import nrPBCHPRBS                  # noqa: E402
from py3gpp.nrPRBS import nrPRBS                          # noqa: E402
from py3gpp.nrSSS import nrSSS                            # noqa: E402

from spectrum_intel.core.seq3gpp import gold, qpsk        # noqa: E402
from spectrum_intel.receivers.nr import ssb               # noqa: E402


def test_gold_and_qpsk_bit_exact():
    rng = np.random.default_rng(0)
    for c in [0, 1, 2 ** 31 - 1, 12345] + list(rng.integers(0, 2 ** 31, 10)):
        n = int(rng.integers(1, 7000))
        assert np.array_equal(np.asarray(_calc_gold(int(c), n)), gold(int(c), n))
        assert np.array_equal(np.asarray(nrPRBS(int(c), 288)), gold(int(c), 288))
    for pci in (0, 3, 17, 1007):
        for ib in range(8):
            assert np.array_equal(np.asarray(nrPBCHDMRS(pci, ib)).ravel(), qpsk(gold(nrPBCHDMRScinit(ib, pci), 288)))


def test_nr_receiver_caches_match_py3gpp():
    for pci in (0, 1, 17, 500, 1007):
        assert np.array_equal(np.asarray(nrSSS(pci), float), ssb._sss_seq(pci))
        for v in range(8):
            assert np.array_equal(np.asarray(nrPBCHPRBS(pci, v, 864)).ravel(), ssb._pbch_prbs(pci, v, 864))
            assert np.array_equal(np.asarray(nrPBCHDMRS(pci, v)).ravel(), ssb._pbch_dmrs(pci, v))


def test_bch_decode_matches_py3gpp():
    rng = np.random.default_rng(1)
    for pci in (0, 17, 1007):
        for lssb in range(3):
            soft = rng.standard_normal(864) * 3
            for x, y in zip(nrBCHDecode(soft, 8, lssb, pci), ssb._bch_decode(soft, lssb, pci)):
                assert np.array_equal(np.asarray(x), np.asarray(y))
