# dji — DJI OcuSync / Lightbridge links and DroneID

**Detects**
- **OcuSync video downlink:** OFDM at 15 or 30 kHz SCS, 10/20/40 MHz, identified by its Zadoff-Chu root.
- **RC uplink:** 72-subcarrier ZC hops, with hop tracking.
- **DroneID:** bursts of 9 symbols at 15.36 MSps (OcuSync 2/3), or 8 (legacy: Mavic 2 and older; the first symbol
  dropped, ZC in symbols 3 / 5).
  - **Frame variant** (research, 2026-10-22, PLAN D4; synthetic only): the CP score is ambiguous by one symbol, so
    the ZC pair is tested at the CP peak and one symbol either side. A CP-correlated OFDM symbol in the 9-symbol
    layout's first slot means 9 symbols; none means 8. Field `droneid_frame`; 8-symbol frames are not
    payload-decoded.
  - **Carrier fix:** the integer bin now comes from the DC dip summed over all the frame's symbols, re-sliced up to
    3 times on the carrier. Bursts between the 3 MHz search slices were up to 35 bins (0.5 MHz) off.
  - Tests: `tests/test_droneid_frames.py` (10 / 3 dB, both variants).

**Decodes** DroneID payloads (serial, GPS) through the external `droneid-decoder-reference` (`droneid_api`, a
separate repo found next to any parent folder or via `$SPECTRUM_INTEL_DRONEID_SRC`). Video and uplink are
identified at the physical layer only (`partial`).

## How
- **Video:** a 15/30 kHz CP-OFDM fold on the 0.5 ms slot, plus a DJI ZC-root search over the sync symbols.
  - Each centre is tried with both the 10 and 20 MHz profiles. Only a link that passes the DJI ZC test (or one
    with the same profile) hides nearby candidates. Otherwise the 20 MHz hypothesis of a 10 MHz link, which folds
    just as strongly but fails the ZC test, hides the right one.
  - The mode (10/20 MHz) comes from the verified ZC length: 601 → 600 subcarriers, 1201 → 1200. The spectrum-edge
    width is kept as `bw_measured_mhz` only, because overlapping Wi-Fi or a far link cuts it short.
  - The carrier estimate re-centres its slice on itself until it stops moving. The strongest scan centre can be MHz
    off the band, and its slice then cuts the band edge and biases the band position.
  - A link that failed the ZC test hides only another hypothesis of itself (same profile, refined carrier within
    0.5 MHz).
  - **Fallback when the ZC test fails** (15 kHz): some links do not keep one symbol grid across bursts. The link is
    then re-timed burst by burst (symbol starts found one by one, up to 3 missing symbols bridged, a leading ZC
    symbol up to P + 256 samples ahead), and the exact ZC test runs on the individual symbols whose DD coherence
    stands > 8 sigma above the rest (`evidence = zc_any_symbol`). No `clock_ppm` for such links.
  - **Robust burst extent under overlapping signals** (Wi-Fi, other links over part or all of the band):
    - the band position comes from the per-bin *median* over symbols, so an interferer in fewer than half of
      them cannot pull the carrier;
    - each symbol is tested on the full band and on 4 sub-bands; it is on when any of them shows the CP
      structure, so Wi-Fi over 3/4 of the band does not cut the burst;
    - a gap of ≤ 0.6 ms in which every sub-band carries at least the power of the neighbouring link symbols
      (the band is covered, the link cannot be observed) is bridged; a gap where any sub-band drops stays;
    - edges: trimmed and extended at 3 per-symbol sigma; a run must contain a symbol where the full-band test
      alone fires (sub-band evidence only extends or bridges a burst, never creates one);
    - RC uplink hops are detected first; video symbols overlapping an in-band hop are excluded (the uplink is
      TDD with the video and has the same numerology);
    - the DJI evidence and the clock fit use only the symbols the full-band test finds on its own.
  - An exact ZC match more than 0.3 above its reference is accepted without the slot statistic (noise
    ~0.1–0.13, Herelink ≤ 0.2).
  - **O4 evidence is positional:** the repeated reference is a comb-4 pilot pair at two slot positions 8 apart
    (4/12 of the 14-symbol slot on D_20260609-151042, 5/13 on the Mini 4 Pro), like NR's TRS. The test scores
    the best of the 14 rotations of {r, r+8} instead of averaging every pair (about 2x the z on every O4 link).
  - **Frame-aware completion at low SNR:** the full-band per-symbol evidence is folded over the TDD frame
    (2.5 / 5 ms) within the link's active span; in every frame whose burst region is significant (z > 4) the
    strongest contiguous stretch is added, its edges trimmed to symbols on by their own evidence; dips of ≤ 6
    symbols inside such a region are joined. Masked-gap bridging runs again afterwards.
  - **Further symbol grids:** a link can re-time mid-capture (D_20260609-151042: the O4 grid jumps by 316 samples,
    10.3 µs, at ~16.8 ms, and the TDD frame phase with it); the capture-level fold keeps one grid. For every reported
    link, further peaks of its own CP fold (z > 8, more than a CP from the grids already used, up to 2) are
    segmented on their own phase at the same carrier and profile; runs not overlapping the link in time (≥ 20
    symbols) are reported as the same link (`link_id`), their own DJI evidence kept in the log.
  - **O4 bandwidth:** the spectrum-edge walk closes dips < 0.6 MHz inside the band (O4's pilot/null pattern):
    the 18 MHz link of D_20260609-151042 measured 7–11 MHz and was labelled 10 MHz; now 18.0 MHz, `video_o4_20MHz`.
  - **Band and bursts from in-band power** (`core/bursts.py`, 2026-10-13): for a link without a verified ZC length
    (O4) the band is the best contiguous rectangle of the on-vs-off power ratio (D_20260609-151042 ch0: 14.1 →
    18.1 MHz, 5761.5–5779.6 MHz); for every link the bursts of all its grids are extended / joined by the in-band
    power runs they touch (sub-band 25th-percentile rank, interferers over part of the band ignored, a stronger signal
    over the band bridged): ch0's split bursts (6.93–7.68 + 8.25–8.60 ms) and late starts are gone, ch0–3 agree within
    ~50 µs. `band_from` says which measurement gave the band.
  - **"-like" tier:** a link whose evidence z is 5–8 (15 kHz: exact ZC test passed) is reported as
    `DJI OcuSync …-like` at confidence ≤ 0.6 (`evidence_tier = like`), one per band and never next to a
    confirmed link, instead of being dropped.
- **Uplink:** the 15 kHz OFDM hop scan and hop tracker in `core.ofdm15k`.
- **DroneID:** ZC sync, then the burst is handed to `droneid_api.decode_ddc_packet`, which does the turbo decode
  and CRC.

## API (depends on `spectrum_intel.core`, plus the optional DroneID decoder)
```python
from spectrum_intel.receivers.dji import detect, droneid_decode
ems = detect(x, fs, fc, {"decode_droneid": True})   # families dji_ocusync / dji_droneid
```

## Validation
`tests/receivers/test_det_dji.py` gives OVERALL PASS:
- Mini 2 video and Mini 4 Pro O4 video;
- D_20260916 uplink (5 hops on 3 channels);
- the DroneID test vectors decode at every added-noise level tested;
- zero false alarms on noise and on the Herelink capture.

Also checked against the detector GUI's boxes on two far Mini 2 captures (2 km, 5770 MHz,
`detection-spectrogram-guis/data/5770_Ant2_Mini_2_*`):
- D_20260916-171218: 10 MHz video at 5760.49 MHz (ZC root 218), 5755.98–5764.98 MHz against the GUI's
  5756.1–5765.2 MHz, every burst within about 0.1 ms;
- D_20260916-171137: 20 MHz video at 5780.48 MHz (ZC root 754), 5771.48–5789.48 MHz against 5771.4–5789.7 MHz.
- Gap: RC uplink hops that fall inside the video band (5764 MHz on 171218) are not found.

Missed links found and fixed on 2026-10-10 (`tmpdata/`, not in testdata):
- **D_20260609-151042 ch0** (FMComms5, 4 channels, drone 2 km at 5770 MHz; an O4-type 30 kHz link although the
  metadata says Mini 2): ch0 sits ~5 dB below ch3 (bursts ~-5 dB in band); the old O4 test scored 6.35 < 8 and
  the link was dropped. Now z 9.5 (ch1–3: 20.9 / 18.9 / 26.5). 2026-10-12: the link re-times at ~16.8 ms, so
  only the 6 bursts after it were reported; with the further-grid search all 14 are (ch3; the second grid's own
  evidence is 24.5, and 5.0 on ch0, where the first half fragments), at 18.0 MHz (`video_o4_20MHz`).
- **D_20260515-192148** (Mini 2, 500 m, 5770 MHz, labelled boxes): the 20 MHz video was rejected because its
  carrier came out 2.35 MHz high (scan peak at 5779 MHz, slice cut the band edge, Wi-Fi at 5785 MHz) and the
  correct 5769 MHz candidate was hidden by the failed one. Now: 20 MHz video at 5770.58 MHz, ZC root 820 (match
  0.878), 7 bursts against the 8 labelled ones (two labelled bursts 0.4 ms apart come out as one).
- **test_signal-20351019 2484 MHz ch3**: a 10 MHz link at 2466.50 MHz under Wi-Fi channel 11. Its grid moves
  10–20 samples per burst (400 samples over 20 ms), so the capture-level grid kept 138 of about 330 symbols and
  none of its 3 ZC symbols (one per ~5 ms burst, root 116). Now found through the re-timing fallback: 8 bursts
  (the 7 visible ones, one split at a half-symbol jump), exact ZC match 0.68 against 0.08.
- On `test_det_dji.py` every testdata emission is unchanged at 0 dB. The carrier re-centring also keeps the
  tv_2p4_122710 10 MHz video to +30 dB added noise (it was lost from +15 dB). tv_5p8_125115 (Mini 2 at 1 m):
  at +10 dB one labelled video burst (14.9–16.1 ms) is now detected, but about 5 MHz high (5775.5 MHz) — a
  single weak burst lets the ZC test pass at a shifted carrier; at 0 dB its video is still not found.
