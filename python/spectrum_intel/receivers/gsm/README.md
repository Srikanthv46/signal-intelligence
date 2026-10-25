# gsm — GSM BCCH carriers (FCCH / SCH)

**Detects** BCCH carriers by their FCCH tone (+67.708 kHz) on the ARFCN raster inside GSM-band energy regions.
**Decodes** the SCH (BSIC and TDMA frame number) when CRC-10 passes (`success`); otherwise it reports the
FCCH-confirmed carrier, ARFCN, CFO and the receiver error measured from the FCCH (`partial`).

## How
1. Each ARFCN is resampled to 4 samples/symbol with the receiver ppm applied (caller / sweep, or the one the `lte`
   receiver measured). The FCCH tone is searched over the carrier offset (`find_fcch_cfo`): 0.4 ms STFT frames
   (2048-point FFT), the largest fraction of frame power in one bin within +67.708 kHz ± 5 ppm (± 30 ppm when no ppm
   is known); `find_fcch` at that tone gives the burst start (95%-plateau) and purity (>= 0.5); the CFO is the
   interpolated FFT peak of the 0.4 ms window from the burst start. (`find_fcch`'s own lag-1 phase CFO is biased by
   kHz by residual adjacent-carrier power in the 220 kHz channel.)
2. A GSM BTS clock is within 0.05 ppm, so every FCCH of one block shows the same receiver error: the median of the
   clean ones (purity >= 0.8, else all) is the reference and a tone more than 1.5 ppm off it is not an FCCH
   (`receiver_ppm_fcch` per carrier). Live (LIVE_940MHz): 13–14 BCCH carriers within 9.56–9.72 ppm, with or without
   the caller's ppm; noise-only channels peak at purity 0.32–0.39.
3. One TDMA frame later: the SCH is demodulated as differential GMSK; the training sequence must match at least
   56/64 within ±20 symbols. Then Viterbi (K=5), CRC-10, and the fields (LSB-first octets: BSIC, T1, T2, T3').
4. Most SCHs still fail the CRC, because of multipath without an MLSE equaliser. That is the next improvement.

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.receivers.gsm import find_fcch, find_fcch_cfo, demod_sch, viterbi, crc10_ok, sch_fields, OFS
```

## Validation
`tests/test_core.py`: the SCH coding round trip (conv. code + CRC-10 + Viterbi with 4 channel errors); an FCCH at
+9.15 kHz (+9.7 ppm) that `find_fcch` misses and `find_fcch_cfo` finds within 30 Hz, and nothing on noise.
Live: `testdata/live/LIVE_940MHz` (GSM-900, 13–14 carriers, BSIC 47 on ARFCN 102 when the SCH passes) and
`testdata/cellular/SWEEP_..._928.08MHz`.
