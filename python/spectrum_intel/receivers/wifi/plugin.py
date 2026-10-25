"""Pipeline plugins for Wi-Fi: OFDM (802.11a/g/n/ac/ax) first, then 802.11b DSSS/CCK."""
from .plugin_dsss import PLUGINS as _DSSS
from .plugin_ofdm import PLUGINS as _OFDM

PLUGINS = _OFDM + _DSSS
