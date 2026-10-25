"""Mode S / ADS-B 1090 MHz: preamble correlation, PPM, CRC-24 (DF17/18 exact, DF11 with
interrogator code), every field of every downlink format, and a per-aircraft summary.

    demod(env, thr) -> list[msg]     env = |x| at modes.FS around 1090 MHz
    fields(msgs) -> msgs             per-frame fields (messages.py, pyModeS)
    tracks(msgs, ref=None) -> list   per-aircraft summary; msgs need 't' seconds (tracks.py)

Depends on spectrum_intel.core and pyModeS. See README.md.
"""
from .messages import fields  # noqa: F401
from .modes import F_ADSB, FS, demod  # noqa: F401
from .tracks import tracks  # noqa: F401
