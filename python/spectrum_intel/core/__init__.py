"""Shared building blocks. Everything here depends only on numpy/scipy; receivers
depend on core and nothing else in the package.

    dsp          pipeline channelizer (absolute-frequency slices), burst masks
    spectral     research channelizer (+ spur notches), spur finder, CP metrics, PSD, spectrogram
    ofdm15k      15 kHz OFDM burst scan + hop tracker
    emission     Emission record + SigMF loaders used by the research receivers
    common       one-line import of emission + spectral for the research receivers
    seq3gpp      LTE PSS/SSS sequences
    bands_india  NFAP-2022 India allocation table, channel rasters, FPV video ranges
"""
