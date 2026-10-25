"""Digital land-mobile radio: TETRA (BSCH: MCC, MNC, colour code), DMR (sync, colour code, data type, BPTC, CSBK
CRC) and NXDN (frame sync), from any IQ block that holds the carrier.

    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.receivers.pmr_digital import detect
    emissions = detect(Channelizer(x, fs, fc))           # list[core.emission.Emission]

Lower level: tetra.demod / tetra.decode_bsch, fsk4.discriminator / fsk4.dmr / fsk4.nxdn.
"""
from .detector import detect  # noqa: F401
