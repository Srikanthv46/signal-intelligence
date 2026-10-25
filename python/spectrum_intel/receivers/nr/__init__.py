"""5G NR SSB receiver: GSCN raster search, PSS, differential SSS, DMRS i_SSB, PBCH
polar decode + CRC-24C (py3gpp) -> MIB.

    _SSBSlice(ch, f_ss, scs), decode_ssb_multi(slice, peaks) -> dict

Depends on spectrum_intel.core and py3gpp. See README.md.
"""
from .ssb import NR_BANDS, decode_ssb_multi  # noqa: F401
