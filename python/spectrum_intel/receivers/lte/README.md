# lte — LTE (E-UTRA) downlink cell search

**Detects** LTE downlink carriers and **decodes** the physical cell identity:
- PSS N_ID2 and SSS N_ID1, giving the **PCI**;
- duplex (FDD/TDD) and cyclic prefix;
- EARFCN and band;
- channel bandwidth from the occupied energy: the largest 3GPP bandwidth whose band around the carrier is
  always-on (>= 85 % of bins > 3 dB) and flat (10–90th percentile spread <= 8 dB), so adjacent always-on carriers
  (GSM beside a 5 MHz LTE) do not widen it;
- the receiver's frequency error in ppm.

Status is `partial`: MIB/SIB (PBCH/PDSCH) are not decoded yet, so there is no operator (PLMN).

## How
1. Candidates: every energy detection ≥ 1 MHz wide and ≥ 3 dB in an LTE downlink band. Continuous ones ≥ 6 dB are
   searched as before; bursty or weaker ones too (a lightly loaded or TDD downlink is bursty in energy, and the PSS
   correlates over many symbols: B28 PCI 351 has 3.6 dB energy SNR and PSS z 77), widest first, skipping rasters
   already searched, and **only an SSS-resolved cell counts** from them. A cell must sit on an LTE band raster
   (Herelink's LTE PHY at 2.41 GHz is left to `digital_vtx`).
2. Carrier hypotheses on the 100 kHz EARFCN raster, corrected by the receiver ppm (caller / sweep, or
   `wide_ppm_search` from the strongest PSS).
3. PSS: a 3.84 MSps slice, 3 roots; offsets 0 Hz when the receiver error is known and its residual (2.5 ppm) is
   < 3 kHz, else ±5 kHz; normalised correlation folded over 5 ms; a hit needs z ≥ 15. A PSS-only hit (old
   candidates only) needs z ≥ 40 and isolation ≥ 3. ZC time-frequency ambiguity at integer-subcarrier offsets is
   removed by NMS.
4. SSS: 168 N_ID1 × {SF0, SF5} hypotheses at FDD/TDD × normal/extended CP, accumulated over half-frames with
   half-frame parity. Each half-frame is **maximum-ratio combined** with the PSS channel estimate smoothed over 15
   subcarriers (225 kHz) on each side of DC (a per-subcarrier zero-forcing equaliser amplified the noise: B5 PCI
   169 went from SSS z 4.1 / ratio 1.43 to 6.9 / 1.42, B28 PCI 351 from undecided to 10.0 / 2.02). It is decisive
   at z ≥ 4.5 and z/z₂ ≥ 1.35.

**Speed** (results identical): the correlation FFT is `next_fast_len(L + N)` long (was the next power of two), the
template spectra are cached per length, the window energy is a cumulative sum, and `core.dsp.Channelizer.slice`
finds its bins by index arithmetic. On LIVE_868 the LTE stage went from 15 s (bursty candidates, before these) to
1.6 s.

## Validation (2026-10-10)
`testdata/live`: LIVE_868 B5 PCI 169 (EARFCN 2450), LIVE_760 B28 PCI 351 (EARFCN 9235) — both missed before;
`testdata/cellular` 770.88 MHz: B28 PCI 415 now found from the single block (before only with the sweep's noise
template); D_20260211-161016 (2381 MHz): B40 TDD PCI 342 (SSS 6.9 vs 3.6). No new PSS-only cells anywhere.

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.core.dsp import Channelizer, IQContext
from spectrum_intel.receivers.lte import wide_ppm_search, search_cells
cal = wide_ppm_search(Channelizer(x, fs, fc), centres_hz)
cells = search_cells(IQContext(x, fs, fc), region, cal["ppm"])   # region: object with f0, f1, fc (Hz)
```
The sequences come from `core.seq3gpp`, which is shared with Herelink.

## Pipeline plugin
Capture-level and owns its band: energy detections that are at least 80% inside the cell's band belong to it.
