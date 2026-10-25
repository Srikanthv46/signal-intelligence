"""GSM BCCH carrier receiver: FCCH tone, SCH (differential GMSK, Viterbi, CRC-10) -> BSIC, FN.

    find_fcch(y) / find_fcch_cfo(y, max_cfo) -> (start, cfo, score); demod_sch(y, start) -> soft bits; viterbi(soft);
    crc10_ok(bits25, par10); sch_fields(bits25)          y at sch.OFS (= 4 x 270.833 kSps)

Depends on spectrum_intel.core only. See README.md.
"""
from .sch import GSM_BANDS, OFS, crc10_ok, demod_sch, find_fcch, find_fcch_cfo, sch_fields, viterbi  # noqa: F401
