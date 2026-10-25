"""DVB-T-framed COFDM links (drone video downlinks), research receiver (research/drone_transmitters/06_dvbt_cofdm.md).

    detect(x, fs, fc, ctx=None) -> list[Emission]      family "cofdm", label dvbt | dvbt_like
"""
from .detector import detect  # noqa: F401
