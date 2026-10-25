# Walksnail header symbols (bench, PLAN D2, 2026-10-22)

Capture: `testdata/live/LIVE_OFDM_VTX_5805MHz` (34 ms, 4 complete frames at 7.163 ms). Slice 22.4 MSps (N 2048,
CP 256) at 5805.07 MHz; frames anchored on the 7.1629 ms grid of the two frames whose CP-located starts agree
(`probe.py`); frame-to-frame subcarrier ratios with the timing ramp and the common (sign-blind) phase removed.

| Symbols (relative to the anchor) | near-real fraction (noise 0.25, clean BPSK ~1) | sign flips frame to frame |
|---|---|---|
| -6 … -4, +3, +7, +8 | 0.22–0.28 | ~0.50 (random: data) |
| -3 … +2, +4 … +6 | 0.30–0.46 | 0.25–0.66 (not constant) |

Reading: a block of ~9 symbols around each frame start carries partly structured, frame-varying content, but not the
clean per-antenna BPSK with a few counter-bit flips the published description implies. Four frames cannot fix a 6-bit
counter's bit positions. **Needed:** a continuous capture of >= 0.5 s (a full 64-frame counter cycle) of a Walksnail
link, ideally with a known 1-antenna and 2-antenna unit. Not added to the library.
