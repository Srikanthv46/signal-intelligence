# emitters — burst → radio → drone-system association

Given Emission records from any receiver (DJI, Herelink, analog VTX, RC links), this decides:
- which bursts come from **one radio**, for example one hopping link seen on several channels;
- how strongly the physics says two radios of the same type are **different**;
- which radios belong to one **drone system**.

## How
- **Timing:** burst starts on one frame-clock lattice (frame period, slot phase; residual rms in ns), from
  CP-train fits re-timed on learned reference symbols.
- **Oscillator fingerprint:** carrier offset in ppm, consistent across hops because one crystal drives them all,
  and clock ppm.
- **Scoring:** a pairwise same-emitter log-likelihood ratio from per-family models (`families.py`), followed by
  average-linkage clustering. Each emitter carries its hop pattern, fingerprint, confidence and evidence.

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.receivers.emitters import measure_fingerprints, associate, _assoc_models
measure_fingerprints(x, fs, fc, emissions)
for em in associate(emissions, fc, models=_assoc_models()):
    print(em.emitter_id, em.family, em.frequencies, em.confidence)
```
`extract_bursts(x, fs, fc)` produces 15 kHz-OFDM bursts directly, for use without the other receivers.

## Pipeline plugin
- `postprocess` adds `report["emitters"]` and tags each signal with `decoded.fields.emitter_ids`.
- `finalize_report` adds the `signal_ids` of each emitter.

## Validation
`tests/receivers/test_drone_assoc.py` (pytest): the 5 DJI uplink hops in D_20260916 form one emitter; injected and
synthetic two-radio cases split correctly; same-clock clones are reported as ambiguous.
