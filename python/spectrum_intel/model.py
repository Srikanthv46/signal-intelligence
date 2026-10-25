"""Records passed between pipeline stages."""
from __future__ import annotations

from dataclasses import dataclass, field

# flags that make a signal a receiver artifact, not an emission (response category "artifact")
ARTIFACT_FLAGS = ("lo_leakage", "iq_image", "receiver_spur", "lo_harmonic_response")


@dataclass
class Box:
    """One connected time-frequency region above the noise floor (absolute units)."""
    t0: float               # s from capture start
    t1: float
    f0: float               # Hz, absolute RF
    f1: float
    peak_snr_db: float
    mean_snr_db: float      # power-weighted mean SNR inside the region
    res: str                # "coarse" | "fine" — which spectrogram found it
    n_cells: int = 0
    kind: str = "burst"     # "burst" | "continuous" (on for >= 80% of the capture)

    @property
    def fc(self) -> float:
        return 0.5 * (self.f0 + self.f1)

    @property
    def bw(self) -> float:
        return self.f1 - self.f0

    @property
    def dur(self) -> float:
        return self.t1 - self.t0


@dataclass
class Signal:
    """One emitter: a set of bursts sharing a frequency span."""
    bursts: list[Box]
    f0: float
    f1: float
    snr_db: float = 0.0
    label: str = "unknown"
    confidence: float = 0.0
    modulation: str = "unknown"
    evidence: str = ""
    features: dict = field(default_factory=dict)
    detection: dict = field(default_factory=dict)   # provenance: how it was found
    decoded: dict = field(default_factory=lambda: {
        "status": "not_possible", "reason": "no decoder run", "method": None, "steps": [], "fields": {}})
    flags: list[str] = field(default_factory=list)   # artifact / quality flags
    classify_trace: list[str] = field(default_factory=list)

    @property
    def t0(self) -> float:
        return min(b.t0 for b in self.bursts)

    @property
    def t1(self) -> float:
        return max(b.t1 for b in self.bursts)

    @property
    def fc(self) -> float:
        return 0.5 * (self.f0 + self.f1)

    @property
    def bw(self) -> float:
        return self.f1 - self.f0
