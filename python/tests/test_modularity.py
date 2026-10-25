"""Loose coupling is a tested property: a receiver's public API must load nothing
from the pipeline or from other receivers, and the pipeline must run with any
receiver switched off."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

from spectrum_intel.decode import registry

ROOT = Path(__file__).resolve().parents[2]                              # <repo>
PIPELINE = ("analyze", "cli", "classify", "decode", "detect", "features", "io", "model", "report",
            "research", "summary", "sweep", "viz")
ANALOG = ROOT / "testdata/Analog_20260804-161815/Analog_20260804-161815.sigmf-meta"


def _loaded_after_import(module: str) -> list[str]:
    code = (f"import sys, json; import {module}; "
            "print(json.dumps(sorted(m for m in sys.modules if m.startswith('spectrum_intel'))))")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True, cwd=ROOT)
    return json.loads(out.stdout.strip().splitlines()[-1])


@pytest.mark.parametrize("name", registry.receiver_names())
def test_receiver_api_is_standalone(name):
    loaded = _loaded_after_import(f"spectrum_intel.receivers.{name}")
    bad = [m for m in loaded
           if any(m == f"spectrum_intel.{p}" or m.startswith(f"spectrum_intel.{p}.") for p in PIPELINE)
           or (m.startswith("spectrum_intel.receivers.") and not m.startswith(f"spectrum_intel.receivers.{name}"))]
    assert not bad, f"{name} drags in {bad}"


@pytest.mark.parametrize("name", registry.receiver_names())
def test_every_receiver_has_a_plugin_and_readme(name):
    assert [p for p in registry.plugins(enable=[name])], name
    assert (ROOT / "python/spectrum_intel/receivers" / name / "README.md").is_file(), name


def test_core_is_leaf():
    loaded = _loaded_after_import("spectrum_intel.core.common, spectrum_intel.core.dsp, spectrum_intel.core.ofdm15k, "
                                  "spectrum_intel.core.seq3gpp, spectrum_intel.core.bands_india")
    assert all(m.startswith("spectrum_intel.core") or m == "spectrum_intel" for m in loaded), loaded


def test_registry_rejects_unknown_names():
    with pytest.raises(ValueError):
        registry.plugins(disable=["no_such_receiver"])


def test_pipeline_with_and_without_a_receiver():
    from spectrum_intel.analyze import analyze
    only = analyze(ANALOG, None, registry.plugins(enable=["analog_vtx"]))
    assert [s.label for s in only.signals if "FPV" in s.label] == ["Analog FPV video (NTSC)"]
    without = analyze(ANALOG, None, registry.plugins(enable=["rc_link"]))
    assert not [s for s in without.signals if "FPV" in s.label]
