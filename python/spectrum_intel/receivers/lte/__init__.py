"""LTE downlink cell search: PSS (5 ms fold) + SSS -> PCI, duplex, CP; receiver ppm
from the PSS carrier vs the 100 kHz raster.

    wide_ppm_search(ch, centres) -> dict                ch = core.dsp.Channelizer
    search_cells(ctx, sig, ppm) -> list[dict]           ctx = core.dsp.IQContext (sig: object with f0/f1/fc)
    sss_decode(y, nid2, instances), sss_decisive(s)

Depends on spectrum_intel.core only. See README.md.
"""
from .cellsearch import LTE_BANDS, search_cells, sss_decode, sss_decisive, wide_ppm_search  # noqa: F401
