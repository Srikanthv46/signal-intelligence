# nr — 5G NR SS/PBCH block and MIB

**Detects** NR SSBs on the GSCN raster:
- n78 at 30 kHz SCS;
- n28/n3/n1/n5/n8 at 15 kHz;
- n40/n41 at 30 kHz.

**Decodes**
- the **PCI**, from PSS + differential SSS;
- the **SSB index**, from the PBCH DMRS;
- the **MIB**, via PBCH polar decode with CRC-24C: SFN, SCS common, k_SSB, CORESET#0 / SearchSpace#0, cell barred.

Status is `success` when the CRC passes.

## How
- Band-wide GSCN search within ±0.49·fs, with offsets of ±12 kHz refined to ±2 kHz. A hit needs
  z ≥ max(40, 2.5 × median).
- SSS is differential (channel-free) because synchronised co-channel cells make PSS-based channel estimates
  useless. Scores are summed over the SSB instances of a beam sweep.
- PBCH:
  - per-SSB DMRS gives i_SSB and a channel estimate;
  - MRC LLRs, descrambled per i_SSB with `nrPBCHPRBS` and summed per half-frame;
  - `py3gpp.nrBCHDecode`.

## API (depends on `spectrum_intel.core` + py3gpp)
`ssb._SSBSlice(ch, f_ss, scs)` and `decode_ssb_multi(slice, peaks)`; see `plugin.py` for the search loop.

## Pipeline plugin
Capture-level and owns its band. The signal span is the union of the host energy and the SSB.
