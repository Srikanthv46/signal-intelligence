"""Plugin discovery. Each sub-package of spectrum_intel.receivers that has a
plugin.py with a PLUGINS list contributes its plugins; nothing is imported until
plugins() is called, and a receiver whose plugin fails to import is skipped (and
reported by unavailable()) instead of breaking the pipeline.

Remove a capability by deleting its receivers/<name>/ folder, or switch it off per
run with ``spectrum-intel run ... --disable analog_vtx`` (names are the folder names).
"""
from __future__ import annotations

import importlib
import pkgutil

from .base import Decoder

# registry order = priority for energy-signal decoding; unlisted receivers follow alphabetically
ORDER = ["wifi", "rc_link", "dji", "analog_vtx", "voice_fm", "pmr_digital", "digital_vtx", "lte", "nr", "gsm", "adsb", "emitters"]

_FAILED: dict[str, str] = {}


def receiver_names() -> list[str]:
    from .. import receivers
    found = sorted(m.name for m in pkgutil.iter_modules(receivers.__path__) if m.ispkg)
    return [n for n in ORDER if n in found] + [n for n in found if n not in ORDER]


def plugins(enable: list[str] | None = None, disable: list[str] | None = None) -> list[Decoder]:
    """Instantiated plugins of the selected receivers, in priority order."""
    names = receiver_names()
    unknown = [n for n in (enable or []) + (disable or []) if n not in names]
    if unknown:
        raise ValueError(f"unknown receiver(s) {unknown}; available: {names}")
    out = []
    for n in names:
        if (enable and n not in enable) or (disable and n in disable):
            continue
        try:
            mod = importlib.import_module(f"spectrum_intel.receivers.{n}.plugin")
        except ModuleNotFoundError as e:
            if e.name == f"spectrum_intel.receivers.{n}.plugin":
                continue                                    # a receiver without a pipeline plugin
            _FAILED[n] = f"{type(e).__name__}: {e}"
            continue
        except Exception as e:                              # noqa: BLE001
            _FAILED[n] = f"{type(e).__name__}: {e}"
            continue
        out += list(getattr(mod, "PLUGINS", []))
    return out


def unavailable() -> dict[str, str]:
    return dict(_FAILED)
