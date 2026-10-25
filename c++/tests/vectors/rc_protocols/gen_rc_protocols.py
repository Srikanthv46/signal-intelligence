#!/usr/bin/env python3
"""Reference for the RC / telemetry protocol research catalog (python rc_link/protocols.py); no capture.
    .venv/bin/python -I c++/tests/vectors/rc_protocols/gen_rc_protocols.py"""
import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link.protocols import candidates   # noqa: E402

cases = [(70.1e3, 56e3, 2440.3, 9.0), (64.2e3, 62e3, 433.9, None), (85.0e3, 41e3, 868.9, 6.67), (250e3, 63e3, 2410, 22.5),
         (31e3, 32e3, 2450, None), (100e3, 50e3, 866.0, None), (1.0e6, None, 2440.0, 2.0), (12e3, 5e3, 433.5, None)]
(OUT / "ref.json").write_text(json.dumps([dict(args=list(a), out=candidates(*a)) for a in cases], indent=1))
