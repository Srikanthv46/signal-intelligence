"""Paths for the receiver validation scripts (moved from wifi-decoding/test_*.py)."""
import os
from pathlib import Path

from spectrum_intel.receivers.dji._paths import droneid_src

ROOT = Path(__file__).resolve().parents[3]                              # <repo>
DATA = str(ROOT / "testdata")                       # D_20260916-171343/, 
OUT = str(ROOT / "out" / "phy_tests")
DRONEID_TV = os.path.join(os.path.dirname(droneid_src()), "test_vectors") if droneid_src() else ""
