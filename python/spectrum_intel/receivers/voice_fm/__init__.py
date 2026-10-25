"""Narrowband FM voice receiver (walkie-talkies, PMR446 / LPD433 / FRS, land-mobile).

    detect(x, fs, fc, ctx=None) -> list[Emission]      family "voice_fm", label nbfm_voice[_<service>_ch<n>]
    VoiceDemodulator.from_fields(fields).process(iq, fs, fc) -> AudioChunk   streaming audio, one chunk per block
    demodulate(iq, fs, fc, fields) -> AudioChunk                          whole capture -> audio
    subaudio: CTCSS tone and DCS code (Golay-verified) from the discriminator output
    synth: synthetic NBFM transmitters (tests)

Analog FM only: scrambled (voice inversion) or digital (DMR / dPMR) voice is not turned
into audio. Depends on spectrum_intel.core only. See README.md.
"""
from .audio import AudioChunk, VoiceDemodulator, demodulate  # noqa: F401
from .channels import CTCSS, DCS_CODES, PLANS, match_channel  # noqa: F401
from .detector import VOICE_RANGES, detect  # noqa: F401
