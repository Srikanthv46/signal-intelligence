"""The shared contract files in schema/ (also used by the C++ implementation):

    response.schema.json   the JSON Schema of a response
    technologies.json      label / receiver -> stable technology code and family
    requirements.json      sample rate and duration each receiver needs
    descriptions.json      plain-language description of every technology (label variants first)

They are read once and cached. ``schema`` next to this module is a link to the
repository's top-level schema/ folder, so both implementations read the same files.
"""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from ..model import ARTIFACT_FLAGS  # noqa: F401  (re-exported)

SCHEMA_VERSION = "1.4.0"


def schema_dir() -> Path:
    here = Path(__file__).resolve().parent / "schema"
    if (here / "technologies.json").exists():
        return here
    for p in Path(__file__).resolve().parents:            # a copy of the package without the link
        if (p / "schema" / "technologies.json").exists():
            return p / "schema"
    raise FileNotFoundError("spectrum_intel schema/ folder not found")


@lru_cache(maxsize=None)
def load(name: str) -> dict:
    return json.loads((schema_dir() / name).read_text())


@lru_cache(maxsize=None)
def _rules() -> list[tuple[str, str, str | None, re.Pattern | None]]:
    return [(r["technology"], r["family"], r["decoder"], re.compile(r["label"]) if r["label"] else None)
            for r in load("technologies.json")["rules"]]


def technology(label: str, decoder: str | None) -> tuple[str, str]:
    """(technology code, family) of a signal from its label and the receiver that defined it."""
    for tech, fam, dec, rx in _rules():
        if dec is not None and dec != decoder:
            continue
        if rx is not None and not rx.search(label or ""):
            continue
        return tech, fam
    return "unknown", "unknown"


@lru_cache(maxsize=None)
def _descriptions() -> list[tuple[str, re.Pattern | None, str]]:
    return [(e["technology"], re.compile(e["label"]) if e.get("label") else None, e["text"])
            for e in load("descriptions.json")["entries"]]


def description(label: str, decoder: str | None, flags=()) -> str:
    """Plain-language description of a signal (schema/descriptions.json): what it is and what the device does. A signal flagged as a receiver artifact is described as that artifact (as
    api.response re-labels its technology). Falls back to the technology family's description."""
    tech, fam = technology(label, decoder)
    if fam != "artifact":
        hit = next((x for x in ARTIFACT_FLAGS if x in (flags or ())), None)
        if hit:
            tech, fam = hit, "artifact"
    for t, rx, text in _descriptions():
        if t == tech and (rx is None or rx.search(label or "")):
            return text
    f = load("technologies.json")["families"].get(fam, "")
    return (f[:1].upper() + f[1:] + ".") if f else ""


def validation(tech: str, decoder: str | None, fields: dict) -> dict:
    """Validation level of a signal (schema/capabilities.json): 'proven' (tested against a real transmitter /
    capture) or 'research' (published specifications, synthetic validation only), plus the decoded fields that
    come from research-level field capabilities. A technology missing from the registry is 'research'."""
    reg = load("capabilities.json")
    t = reg["technologies"].get(tech)
    out = {"level": t["level"] if t else "research"}
    if t and t.get("proof"):
        out["proof"] = t["proof"]
    if t and t.get("sources"):
        out["sources"] = t["sources"]
    rf, src = [], []
    for c in reg["field_capabilities"]:
        if c["level"] != "research" or c.get("decoder") not in (None, decoder):
            continue
        hit = [k for k in c["fields"] if fields.get(k) not in (None, [], {}, "")]
        if hit:
            rf += [k for k in hit if k not in rf]
            if c["sources"] not in src:
                src.append(c["sources"])
    if rf:
        out["research_fields"] = rf
        out["research_sources"] = src
    return out


def identifier_keys() -> list[str]:
    return load("technologies.json")["identifier_keys"]


def requirements() -> dict:
    """Per-receiver sample-rate / duration needs, and the pipeline's minimum block size."""
    return load("requirements.json")
