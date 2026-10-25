# Live LO-by-LO bench sweep (Pluto, usb:1.2.5)

Setup: 61.44 MSps, rf bandwidth 56 MHz, 2^21 samples per buffer, manual gain 40 dB (stepped down 10 dB if
>0.1 % of samples clip), receiver_ppm 7.4. LO step 36 MHz (56 MHz analog band, 20 MHz overlap), first LO 328 MHz
(300–356 MHz). Per LO: rx buffer destroyed before the retune, 4 buffers discarded, 3 kept 0.3 s apart.
Tools: `lo_step.py` (capture + Engine), `inspect_lo.py` (spectrogram + zooms), `bursts.py` (burst measurements).

## Capture-chain findings
- **Stale buffer after retune:** buffer 0 differs (median PSD change 2.4 dB into buffer 1); from buffer 2 the
  spectrum is stable (0.14–0.22 dB). Discarding 4 is enough.
- **Self-noise step inside every buffer:** at ~12 ms (sample ~740k of 2^21) the floor rises 4–5 dB across the
  whole 56 MHz and stays up to the end of the buffer; same sample offset in all three snapshots, so it is the
  capture chain (USB transfer / host), not the air. Also the first ~0.5 ms is 5 dB quieter. Costs ~4 dB of
  sensitivity on two thirds of each buffer near 300–360 MHz. Acquisition-side fix (cable/ferrite/antenna
  away from the Pluto and the laptop) - to check at other LOs.
- **Self spurs:** 307.2 MHz = 5 x 61.44 MHz (vanishes at 50 MSps); 320.0 MHz stays at 50 MSps (8 x 40 MHz
  reference likely). The Engine labels both "PMR land-mobile (narrowband)-like" (conf 0.55) -> false label.

- **LO-harmonic images (the big one).** The Pluto has no RF preselector; its mixer also converts through odd LO
  harmonics. A strong signal at F appears at baseband d = 3LO - F (H3, inverted) or d = F - 5LO (H5) etc.
  Test (`harmonic_map.py`, automatic in `lo_step.py --dither 1`): capture LO and LO+delta; real signals stay
  put, H3 images move +4 delta, H5 -4 delta, H7 +8 delta. At this site the 1.8 GHz cellular band (LTE B3 20 MHz
  carrier ~1817.5-1835.5 MHz + GSM-1800 hopping bursts) fills 336-392 MHz via H5, and GSM-900 BCCH/hopping
  (955-959 MHz) shows up at 351-357 MHz via H3. The Engine reports them as `unknown` / "TETRA-like" /
  "PMR-like" at the wrong frequency. Hopping bursts are matched poorly by mean spectra; steady carriers are
  the reliable markers. Hardware fix: a low-pass/band-pass preselector per sub-band below ~1.2 GHz.

## Per LO
### 328 MHz (300–356)
- (RESOLVED: images, see below) 352.6–355.0 MHz "frequency hopper": bursts 0.550 ms exactly, ~200–240 kHz wide, 2-FSK/GFSK (instantaneous
  frequency ±70 kHz), a hop every 1.7–2.9 ms. Engine: many separate `unknown` detections, no emitter grouping.
- tones 355.86 and 356.6 MHz (+22 / +18 dB) at the band edge; 352.45/352.98 CW-like weak.
- nothing decoded.
- **Resolution:** tones 355.86/356.6 MHz move -4 MHz when LO goes 328->327: H3 images of GSM-900 carriers at
  956.14/955.4 MHz. The "hopper" (0.55 ms, ~200 kHz, ±70 kHz FSK = GSM GMSK bursts) is GSM-900 hopping traffic via
  H3 (957-959 MHz). Fixed in absolute frequency (not images): 300.0 MHz CW (+20 dB, origin unknown: 12 x 25 MHz
  clock?), 307.2 MHz (5 x 61.44 MHz sample clock), 320.0 MHz (8 x 40 MHz reference likely).
- Real on-air signals 300-356 MHz at this site: none found besides those CW lines.

### 364 MHz (336–392)
- Everything structured is an H5 image (LO 364 vs 365: every segment shifts -4 MHz, corr 0.85-0.99): LTE-like
  15 kHz OFDM block 362.5-379.5 MHz (CP lag 66.67 us, 1 kHz subframe lines; PSS smeared, no decisive SSS) =
  LTE B3 carrier at ~1817.5-1835.5 MHz; GMSK bursts at 358.5-361.5 / 379.5-381.5 = GSM-1800 at 1815-1817 /
  1836 MHz. Engine "TETRA-like" at 391.5/392.0 and "PMR-like" at 387.7 are on top of image content.
- Self-noise step at ~12 ms again.

### 400 MHz (372–428)
- Harmonic map: 389–393 and 417–423 MHz stay put (not images). But the **416–424 MHz comb** (~0.43 MHz spaced lines,
  420.45 strongest) exists only in the first ~12 ms of every buffer (14 dB there, 1 dB after) -> **self-generated**,
  locked to the capture (third test needed besides the dither map: buffer-time locking). The broadband rise at
  372–390 MHz after 12 ms is the same self-noise step. Engine "PMR-like" at 416.95/417.66/421.97 sit on that comb.
- **TETRA, decoded.** 390.0297 and 391.527 MHz: ~22 kHz wide, 18,000.2 Hz symbol-rate line, pi/4-DQPSK (4
  differential-phase clusters 90 deg apart). Dwell 1.75 s at 2.4 MSps (390.78 MHz, gain 50, `dwell.py`):
  - 390.03 MHz = **main carrier**: 47 sync bursts (38-bit y training, all at the same slot phase mod 510 bits);
    BSCH decoded (`tetra_bsch.py`, from EN 300 392-2: e=0 scrambling, K=120/a=11 deinterleave, RCPC 2/3 ->
    rate-1/4 K=5 Viterbi soft, CRC-CCITT): **39/47 CRC ok (hard decisions 25/47): MCC 404 (India), MNC 4040,
    colour code 7, system code 2**, TN/FN/MN count consistently (MN15 -> 16). Validation: CRC + consistent
    frame counters + India MCC. No external oracle used (osmo-tetra would be one; AGPL -> oracle only).
  - 391.527 MHz: TETRA modulation, no sync bursts -> secondary/traffic carrier of the same network (likely).
  - Engine: labels them `unknown` / "TETRA-like" (0.55), "not_possible". A TETRA receiver needs ~1.1 s at a
    narrow rate (BSCH every ~multiframe on the MCCH; here SBs appear in many slots) — 34 ms blocks are too short
    for reliable SNR; at 61.44 MSps a 2^21 block holds ~600 symbols = ~2 slots.
- 399.994 MHz LO leakage (Engine correct). 430.08 MHz pure CW (+10 dB, steady) at the edge — next LO.
- Next: TETRA SYSINFO (BNCH on the MCCH: main carrier, band, offset, duplex spacing, LA) needs the cell
  scrambling code from MCC/MNC/CC — possible follow-up.

### 436 MHz (408–464)
- H5 images 409–415 MHz (= LTE B1 DL ~2154–2158 MHz).
- **False LoRa/ExpressLRS:** 25 "LoRa sf6 bw500k" / "ExpressLRS-like" detections at 424–429 MHz, all between
  0.6 and 12.3 ms of the buffer, i.e. inside the capture-locked self-noise (diagonal chirp-like stripes there).
- **Self spur:** 430.08 MHz = 7 x 61.44 MHz (raw 430.079996, no ppm offset) -> Engine "Analog FM voice" is false.
- **EMI comb of a 12 MHz clock (+14.4 ppm):** lines at 408.003, 420.003 (+31 dB), 432.003, 456.003 MHz, each with an
  +8 kHz companion; the source-clock error from each line agrees (14.3/14.4 ppm). 12 MHz = USB full-speed
  crystal -> a USB device near the antenna. Engine "Analog FM voice (420.0029 MHz)" is false (FM: 1 kHz rms
  deviation, no audio).
- CW lines 440.0 (raw 439.99804, i.e. +3 ppm source - another clock, 11 x 40 or 110 x 4 MHz) and 433.335 MHz
  (drifting -90 Hz/s, free-running oscillator).
- **456.05 MHz**: real FM carrier ~8 kHz wide (PMR band), in-channel SNR 9 dB even at gain 70 (floor is external
  noise, not ADC) -> not decodable here.

### 472 MHz (444–500)
- 452–481 MHz: **LTE TDD (B40, 2338–2365 MHz) via H5** - clear 5 ms TDD frame pattern in the spectrogram.
- 481–494 MHz: **capture-locked self noise**: diagonal stripes before 12.3 ms, a ~0.43 MHz line comb after;
  fixed absolute frequency (so the dither test says "real") but the in-band spectrum shape changes at the
  step (r = 0.03–0.22) -> host/USB EMI driven by the capture. `triage.py` now tests power and shape locking.
- 480.0034 = 40 x 12 MHz EMI comb (same +7.0 ppm vs our display).
- Real, weak: 456.05 MHz land-mobile FM (seen 2/3); 446.05 MHz PMR446 ch 5 (446.05625; seen 1/3, SNR 4 dB);
  500.003 MHz weak line (4 MHz grid?).

## Triage tool (`triage.py <LO>`), per Engine detection group
SELF-window (t < 12.4 ms) / SELF-clock (n x 61.44 or n x 40 MHz, raw) / SELF-locked (power or spectral shape
differs before/after the 12.3 ms step) / EMI (site 12 MHz comb) / IMAGE Hn (dither) / REAL.

### 508 MHz (480–536) — IQ-only analysis (`modclass.py`)
- **533.336 MHz (+19 dB) = Pluto DDR3 clock** (16 x 33.333 MHz = 533.33 MHz, raw +4 ppm = Zynq PS oscillator, not
  the AD9363 TCXO): constant envelope (kurtosis 1.01), autocorrelation 0.93 at 896 us (unmodulated), 99 % of
  power within ±60 kHz. Lines at 533.335 ± k x 1.92 MHz (527.62, 529.50, 531.42, 535.26, 537.18) belong to it.
- 479–492 MHz capture-locked self noise again; 528.004 = 44 x 12 MHz EMI.
- 510.75 / 519.98 / 527.6 MHz detections: Gaussian envelope (kurtosis ~2.1–2.4), in-band SNR ~ -3 dB: noise.
- No UHF TV (DVB-T2 / analog PAL) in 480–536 MHz at this site.

### 544 MHz (516–572)
- More own-clock lines: **Zynq PS 33.333 MHz harmonics** (500.0027 = 15x, 533.336 = 16x DDR, 566.670 = 17x; all
  +5 ppm) and lines at exact integer MHz in our display (536.000016, 552.000011: within 16 Hz -> derived from the
  Pluto reference, since an external line would sit ~4 kHz off). Triage covers both.
- **Second EMI comb: 40 MHz clock at +2.9 ppm** (source error from 440.0 / 520.0 / 560.0 MHz lines: 2.95 / 2.89 /
  3.03 ppm). 519.998 MHz = carrier + 14.6 kHz companion, 0.95 kHz rms deviation, no audio -> Engine "Analog FM
  voice (519.9977 MHz)" is false. 560.0 MHz pure CW (kurtosis 1.13).
- The background noise is amplitude-modulated at 7,996 Hz and 96,004 Hz (12 x 8 kHz) in every window:
  **USB 2.0 high-speed microframe rate (125 us)** -> the capture self-noise is the USB link.
- 563.624 MHz: present in the 61.44 MSps snapshots (SNR ~9 dB in 66 kHz, kurtosis 1.5), gone in the 1.75 s dwell
  -> intermittent, not identified.
- No real decodable emission 516–572 MHz.

### 580 MHz (552–608)
- 597.4–602.6 MHz: one comb, lines every ~0.42–0.43 MHz centred on 599.997 MHz = 15 x 40 MHz of the +2.9 ppm EMI
  clock (same 0.43 MHz sideband spacing as the 420 MHz comb) -> EMI, not emissions (triage splits it into
  REAL/SELF-locked pieces; it is one source).
- 552.96 / 566.67 own clocks, 560.0 EMI; 569.58 and 594.78 MHz single weak lines (1/3 snapshots).
- No real decodable emission 552–608 MHz.

### 616 MHz (588–644) and 652 MHz (624–680)
- 617–629 MHz: H3 images of LTE B3 (1835–1846 MHz).
- 585.5 / 586.75 MHz "480 kHz signals": in a 1.75 s dwell the band is flat (±2 dB) and its power is modulated at
  8 / 16 / 24 / 32 kHz (USB microframe) -> USB EMI, not emissions. Same for the 597–603 comb sidebands.
- 640.0 = 16 x 40 MHz EMI (+2.9 ppm clock), 660.004 = 55 x 12 MHz EMI; 658.3/659.6 = its sidebands.
- 620–624 MHz narrow lines: 1 of 3 snapshots, ~6 dB -> too weak/brief to classify.
- No real decodable emission 588–680 MHz.

## Summary 300–680 MHz
Decoded: TETRA 390.03 MHz (MCC 404, MNC 4040, CC 7). Real but weak: 456.05 MHz land-mobile FM, PMR446 ch 5,
CW lines 433.335 / 440 / 560 MHz (EMI clocks). Everything else the Engine reports here is a harmonic image of
1.8/2.1/2.3 GHz cellular, a Pluto clock line (61.44 / 40 / 33.333 MHz families, DDR 533.33), an external clock
comb (12 MHz +14.4 ppm; 40 MHz +2.9 ppm), or USB-link EMI (8 kHz microframe modulation, locked to the capture).

### 688 / 724 MHz (660–752)
- No real emissions: capture-locked noise, own clocks, EMI; 731–743 MHz = H3 images of LTE B1 (2154–2164 MHz).

### 760 / 796 MHz (732–824) — cellular B28 / n28
- **5G NR n28, SSB 782.55 MHz (NR-ARFCN 156510): two cells, PCI 61 and PCI 98** (Engine decoded, MIB CRC ok in
  all 6 snapshots: SCS 15 kHz, ssb_SubcarrierOffset 6, CORESET0 6 / SS0 5, not barred, SIB1 scheduled). SFN
  advances 38 frames per 0.38 s between snapshots (930 -> 968 -> 1006): consistent. At LO 796 PCI 98 is the
  primary with PCI 61 as co-channel; at LO 760 only PCI 61 -> both real, relative level varies.
- **LTE B28, ~760.5 MHz (EARFCN 9235), PCI 351, FDD normal CP, ~5 MHz** (757.3–763.3 MHz real energy):
  bench PSS/SSS (`lte_probe.py`): snap0 decisive (SSS z 6.74 vs 4.0), snap1 same PCI (z 5.9 vs 4.95), snap2 none.
  **Engine missed it** ("LTE-like downlink B28/n28", not_possible): one 34 ms block = 6–7 PSS, SSS marginal.
  Library improvement: accumulate PSS/SSS evidence across blocks of one tuning (StreamAnalyzer / snapshots).
- Edge trap: at LO 796, content beyond ±30.72 MHz (e.g. 758 MHz) aliases (758 -> 819.4 MHz) — the analog filter
  (56 MHz) does not remove it fully; 765.9 MHz (+20 dB) at offset -30.1 MHz is the alias of B5 uplink 827.3 MHz.
  Analyse only within ±28 MHz of the LO.
- 806–824 MHz at LO 796: H3 images of LTE B40 (2361–2370 MHz).

### 832 / 868 MHz (804–896) — B5 and 860 MHz land mobile
- **LTE B5, 874.0 MHz (EARFCN 2450), PCI 169, FDD normal CP, 10 MHz (869–879 MHz)**: bench PSS/SSS decisive in
  2/3 snapshots (SSS z 6.57 vs 3.74; 6.08 vs 3.82), receiver error +9.5–9.7 ppm.
  **Engine missed it: the detector labels the 9.84 MHz carrier `burst` (lightly loaded LTE DL + capture noise
  step) and `lte/plugin.py capture_signals` only searches candidates with a `continuous` burst.** Bench fix
  (`lte_fix_test.py`, in-memory patch: drop the `continuous` requirement, keep the SSS-decisive gate):
  Engine as is 0/3 -> patched 2/3 PCI 169 + 1/3 PSS-only. The B28 PCI 351 miss is very likely the same cause.
  -> candidate library change; must be validated on testdata (false positives) before python/.
- **NXDN, 859.3455 MHz**: 4FSK, 4800.0 sym/s, 12.5 kHz, channel SNR 18.5 dB, NXDN frame sync word (0xCDF59) every
  384 bits = 40 ms frame (42 hits in 1.75 s) (`dmr_probe.py`). Not decoded further.
- **DMR base station, 859.8217 MHz**: 4FSK 4800.1 sym/s, BS-data sync every 288 bits (58 bursts / 1.75 s, 0 sync
  bit errors with the per-burst RRC + sync-aided demod in `dmr_csbk.py`). Slot type = exact Golay(20,8)
  codeword: **colour code 0, data type CSBK** (all bursts). BPTC(196,96) verified from the spec (ETSI TS 102 361-1
  B.1.1, interleave x181): 516/522 rows Hamming(15,11) ok, 860/870 columns Hamming(13,9) ok. Two CSBKs alternate
  (99 00 08 00 e6 40 bb 00 00 00 | f3 9d and 99 10 08 00 e6 40 b9 00 00 00 | e5 f6) but **neither passes the
  standard CRC-CCITT with the CSBK mask A5A5** (CRC routine checked against the XMODEM vector; no constant mask
  fits) -> non-standard/vendor control channel; opcode not claimed.
- 840.006 = 70 x 12 MHz EMI; 806–824 MHz H3 images of B40; 765.9 MHz at LO 796 = alias of B5 uplink.

### 904 MHz (876–932) — GSM-900 / B8 uplink (880–915 MHz in India)
- **GSM handset uplink, 899.6 MHz = ARFCN 48** (DL pair 944.6 MHz): 0.58 ms bursts, 200 kHz, GMSK frequency
  swing +-67.7 kHz; `gsm_ul.py` (GMSK differential demod at 4 sps, normal-burst TSC match): **TSC 3 with 0/25
  errors in all 7 bursts** (2 snapshots; none in the 1st = DTX/idle) -> serving cell BCC 3. Bursts 1 or 2 TDMA
  frames apart (4.615 ms grid). Engine: `unknown` (the GSM receiver handles only the downlink FCCH/SCH).
  Trap: TSC 0 inverted equals TSC 3 shifted one symbol -> match with the fixed (TS 45.004) polarity only.
  -> candidate library feature: GSM uplink burst detector (GMSK, 577 us on the 4.615 ms grid, TSC -> BCC).
- 903.5–906.1 MHz: LTE B8 uplink-like (SC-FDMA bursts 0.1–1 ms, up to 400 kHz) — handset, not decodable.
- 917.7–918.7 MHz: weak narrow bursts (uplink-band, unidentified).
- Lines 880/920 (40 MHz family), 888/900/912 (12 MHz family) = external clock combs; 921.6 = 15 x 61.44 (own).
- 876–879.3 MHz: upper edge of LTE B5 PCI 169 (874 MHz, decoded at LO 868), truncated here.

### 940 MHz (912–968) — B8 / GSM-900 downlink (935–960 MHz), ADC near clipping (0.04–0.06 %) at 40 dB
- **LTE B8 FDD PCI 96, EARFCN 3608 (940.8 MHz), 5 MHz (938.3–943.3)** and **LTE B8 FDD PCI 101, EARFCN 3690
  (949.0 MHz)** (energy 946–953.5 MHz): Engine decoded both in 3/3 snapshots (SSS z 10.6–11.4 vs 5.4–5.6).
  Receiver error from LTE: **+9.5 to +10.1 ppm** (TCXO drifted from the 7.4 ppm calibration).
- **GSM-900 BCCH carriers** (`gsm_dl.py`, library FCCH/SCH on the 200 kHz raster, --ppm 10):
  - **ARFCN 106 (956.2 MHz): SCH CRC OK in snap2 -> BSIC 47 (NCC 4, BCC 7), FN 1777973**; FCCH 3/3 (+9.5 ppm).
  - ARFCN 97 (954.4) and 102 (955.4): FCCH 3/3 (+8.8–9.3 / +6.8–7.0 ppm), SCH training 58–62/64 but CRC fails
    (differential detection without an equaliser — the known MLSE gap). 102 and 106 have FCCH at the same instant
    -> same BTS. ARFCN 11/12 (937.2/937.4): GSM carriers, FCCH weak/odd.
  - ARFCN 45–54 (944–946 MHz): strong 200 kHz carriers; "FCCH" hits there have inconsistent ppm (0–6) and no SCH
    -> not trusted. ARFCN 48 (944.6) is the DL pair of the TSC-3 handset at LO 904 (TCH carrier).
  - **Engine missed every GSM cell here**: it applies the fixed +7.4 ppm; the 2.4 kHz residual rotates the
    FCCH tone by ~1 cycle within the 0.4 ms purity window, so find_fcch fails. With 10 ppm it finds them.
    -> library fix candidate: GSM FCCH search over a CFO grid (or reuse the LTE/NR-measured ppm of the block),
    plus the MLSE equaliser for SCH.
- 958–960 MHz: GSM carriers (lower power); 960+ empty. 935.4/942.3 MHz: IQ images (artifact, correct).

### 976 MHz (948–1004) — top of GSM-900 DL, start of the aeronautical band (960–1215 MHz)
- With `--ppm 10` the Engine now reports GSM BCCH (FCCH) on ARFCN 96–98, 100–102, 106, 108, 109 (partial, SCH not
  decoded). `gsm_dl.py`: real FCCH (ppm +8.6–9.5, purity >= 0.93, SCH training 57–63/64) on ARFCN 97, 100, 106;
  ARFCN 102 consistently +6.7–7.0 ppm (2.6 kHz low: different BTS clock?). **False FCCH**: ARFCN 98 (-13/-14 ppm,
  purity 0.75), 101 (+18 ppm, 0.57), 109 (-0.9 ppm), coinciding in time with the real FCCH of 102/106.
  -> GSM receiver fix candidates (with the 940 ones): FCCH tone search over a CFO grid; accept an FCCH only when its
  residual CFO agrees with the block's ppm (+-2 ppm) and purity >= 0.8 or the SCH training matches; MLSE for SCH.
- 960–1004 MHz: no DME/aeronautical pulses above the noise; 983.04 = 16 x 61.44 (own), 1000.0 integer-MHz line.
  1005–1007 MHz is at the edge (seen at LO 1012).

### 1012 MHz (984–1040) — aeronautical (DME, SSR 1030)
- No emissions: no DME pulse pairs, no 1030 MHz interrogations. Lines: 1000.0 / 1008.0 / 1020.0 MHz (integer-MHz,
  12 MHz family), 1033.338 = 31 x 33.333 MHz (Zynq), 1026.29 weak narrow (unidentified, 2/3 snapshots).
- The strong 1005–1007 MHz lines seen at LO 976 are gone here: they were aliases (f - 61.44 MHz ~ 945 MHz GSM)
  past the 56 MHz analog filter edge (the known edge trap).

### 1048 MHz (1020–1076)
- No emissions. 1066.66 MHz (540 kHz wide, 15 dB, real per dither) = 2 x 533.33 MHz (Zynq DDR clock), with lines
  1.92 MHz (61.44/32) apart: 1063.15, 1064.75, 1068.59, 1070.51 MHz. 1033.34 (31 x 33.333), 1040 (40 MHz family),
  1056 (12 MHz family). 1046–1049.4 MHz moves with the LO (dither) -> receiver artifact around the LO, reported as
  `unknown` 6–11 dB by the Engine (not flagged).

### 1084 MHz (1056–1112) + 1.05 s dwell at 1090 MHz (8 MSps, 4 MHz, gain 60)
- The 34 ms LO-step blocks show no Mode S (expected: sparse squitters; requirements say 10 s recommended).
- **Dwell (`dwell.py --freq 1090 --rate 8e6 --n 8388608`) -> Engine: Mode S / ADS-B success, 17 frames with valid
  parity (DF17 x8, DF11 x4, DF0 x4, DF20 x1), 4 aircraft, all India (0x80xxxx):** 8015CB 11,400–11,425 ft
  (Comm-B DF20), 8015D1 12,550 ft (ADS-B v2), 80175F 32,975 ft level, 80028D 9,500 ft (level-1 transponder,
  DF0/DF11 only). No positions (an even/odd CPR pair needs up to 10 s). The signal's snr_db is reported as 0.0 (pulsed).
- **False positive: "Analog FPV video (PAL)", decode success** at 1089.377 MHz (sync tip 1090.043 MHz = the LO,
  score 19.7, SNR -18 dB): line rate 15625.007 Hz = exactly fs/512 (8 MSps) at +0.46 ppm, while any real
  transmitter is seen ~-10 ppm off through the Pluto's clock error -> a capture-locked artifact at the LO.
  -> analog_vtx fix candidate: reject a line rate locked to fs/N (within ~2 ppm) when receiver_ppm is known/non-zero
  and/or a sync tip at the LO; also SNR < 0 dB with "success" should not happen.
- Lines: 1066.66 (DDR x2) + 1.92 MHz side lines, 1080 (40/12 MHz families), 1100.0, 1104.0, 1105.92 = 18 x 61.44
  (own, flagged), 1107.3 (unidentified narrow).

### 1120 MHz (1092–1148)
- No emissions: 1100.0 / 1104.0 / 1140.0 integer-MHz lines, 1105.92 = 18 x 61.44 (own), weak narrow unidentified
  1128.11 / 1134.31 / 1107.3 MHz (3/3 snapshots, 6–7 dB).

### 1156 MHz (1128–1184) — H3 image of 5G n78, decoded through the image
- Buffer RMS 50–119 vs ~8 at the neighbouring LOs; broadband (whole 56 MHz, 10–21 dB above the floor), on 4 ms /
  off 1 ms every 5 ms (TDD). Not a settling transient (8 discards: same).
- H3 image (library formula: source = 3 LO - (F - LO), inverted) -> 3440–3496 MHz = n78. **Conjugated IQ analysed
  as centred on 3 x LO = 3468 MHz -> Engine: 5G NR n78 PCI 972, SSB 3449.28 MHz (NR-ARFCN 629952), MIB CRC ok 3/3,
  SFN 594 -> 632 -> 668** (consistent with the 0.36–0.38 s between snapshots), SCS 30 kHz, ssb_SubcarrierOffset 12.
  (A first try with a 2 x LO centre gave only "OFDM 30 kHz SCS" -> wrong mapping; the library formula is right.)
- Engine: every signal here carries `possible_lo_harmonic_image` (band-plan warning works). But it also lands on the
  integer-MHz lines 1140/1160 and the 19 x 61.44 = 1167.36 MHz spur, which loses its `receiver_spur` label in the
  snapshot where the image raises the floor (labelled 2/3 only).
  -> library feature candidate: **confirm an image by decoding it** (cellular receivers on conj(x) at 3 LO and on x
  at 5 LO for flagged broadband energy) -> "image of NR n78 PCI 972" instead of a caution; and keep the spur label
  when the floor rises.

### 1192 MHz (1164–1220) — H3 image of a second n78 block (3548–3600 MHz)
- Same 5 ms TDD broadband energy (4–19 dB above the off periods), from 1168 MHz up = H3 sources 3600 MHz down to
  3548 MHz (edge exactly at 3600 MHz -> an n78 carrier ending at 3600). `image_decode.py 1192 --n 3`: no SSB in the
  window (the carrier's SSB lies outside 3548–3604). Nothing else real: 19 x 61.44 spur (labelled 2/3 again).
- Tools: `image_decode.py` (Engine on conj(IQ) at 3 LO / IQ at 5 LO), `floor_profile.py` (per-ms power + on-off
  spectrum), `lo_summary.py` (signals grouped across snapshots).

### 1228 MHz (1200–1256) — GPS L2 / DME top / L-band radar
- No emissions (GPS L2 below the noise as expected). Lines: 1228.80 = 20 x 61.44 (own), 1240.0 / 1252.0
  (reference-locked, flagged), 1242.31 weak. n78 image floor gone (RMS 13).
- **False "DME/TACAN pulse pairs-like" (conf 0.5, 1/3)** at 1200.88 / 1201.33 MHz (band edge): the "pulses" are
  narrow (< 35 kHz) one-STFT-frame blips on a 0.44 MHz comb (1200.44 ... 1202.20). Real DME: 3.5 us pulse pairs
  ~1 MHz wide on the 1 MHz channel raster -> classifier fix candidate (require pulse width / bandwidth / raster).

### 1264 MHz (1236–1292) — L-band radar / amateur 23 cm
- No emissions (no radar pulses in 3 x 34 ms). Lines: 1290.24 = 21 x 61.44 (own), 1240.0 (labelled reference
  spur at LO 1228, `unknown` here), 1269.29 weak.

### 1300 MHz (1272–1328)
- No emissions. 1290.24 = 21 x 61.44: `receiver_spur` 2/3, `unknown` 1/3 (detected 33.8 kHz wide) -> the spur label
  depends on the detected width; 1320.0 integer-MHz line; 1296.3 weak (1/3).

### 1336 MHz (1308–1364)
- No emissions. 1320.0 (integer MHz), 1333.338 = 40 x 33.333 (Zynq), 1351.68 = 22 x 61.44 (spur label 2/3 again).

### 1372 MHz (1344–1400)
- No emissions. 1351.68 spur (3/3 labelled), 1360.0 integer MHz. 1371–1373 moves with the LO (around-LO artifact,
  below the Engine's 6 dB report threshold).

### 1408 MHz (1380–1436) — radio astronomy 1400–1427 (protected)
- No emissions. 1400.0 integer MHz, 1413.12 = 23 x 61.44 (labelled 3/3), 1404.36 weak; around-LO artifact.

### 1444 MHz (1416–1472) — L-band IMT 1427–1518
- **Unidentified real emission 1439.19 MHz** (dither: real), continuous over the block, ~75–90 kHz wide, 8.6 dB, 2/3
  snapshots; envelope Gaussian-like (kurtosis 2.1), |x|^2 line at +-96.0 kHz (= 61.44 MHz / 640, so possibly
  capture-related), weak x^2 line. Engine: `unknown` (correct). 1440.0 integer line; 1433.34 = 43 x 33.333.

### 1480 MHz (1452–1508)
- 1485.42 MHz real (dither), ~32 kHz, 10 dB, 2/3: frequency wobble at 47.3 Hz (+ harmonic 94.6 Hz) in the FM
  baseband, no symbol-rate line -> EMI-like drifting oscillator, unidentified (Engine `unknown`, correct).
  1474.56 = 24 x 61.44 (labelled 3/3), 1500.0 integer line.

### 1516 MHz (1488–1544) — start of MSS downlink (Inmarsat 1525–1559)
- No decodable emissions. 1536.0 = 25 x 61.44 (own, labelled); 1538.568 = unmodulated tone (x^2 line at DC; = 1536 +
  2.56 MHz = 61.44/24 -> clock intermod likely); 1539.49 weak, no structure; 1512.45/1512.50 real (dither) weak
  narrow pair, no symbol line; 1520.0 integer line. No Inmarsat-like 1200 sym/s BPSK found.

### 1552 MHz (1524–1580) — MSS downlink / GNSS L1 (1559–1610)
- No emissions (GNSS below the noise, as expected). 1538.56 / 1539.49 seen again at the same RF -> fixed lines;
  1560.0 integer line; 1561.6 weak (1/3).

### 1588 MHz (1560–1616) — GNSS L1 / GLONASS L1
- GNSS below the noise (as expected). **1599.97 MHz strong real line** (dither peak 28 dB, 15.6 dB Engine SNR,
  ~1.2 MHz skirt, 3/3): tone (x^4 line at DC, lag correlation 0.8) with +-30 kHz FM-baseband lines; ~1600.000 MHz
  after the +10 ppm receiver error -> 40 x 40 MHz (the external 40 MHz comb) with spread-spectrum/FM clocking.
  Engine `unknown` (an inferred "clock harmonic" label would need the comb model).
- Narrow pairs 36 kHz apart recur 54 MHz apart: 1512.45/1512.50 and 1566.46/1566.49 (+ 1538.56 midway?) -> an
  external comb, unidentified. 1576.0, 1560.0 integer-MHz lines.

### 1624 MHz (1596–1652) — Iridium 1616–1626.5
- **External 27.008 MHz comb** (resolves the "pairs" at 1512–1647 MHz): 1512.47 / 1539.48 / 1566.48 / 1620.49 /
  1647.51 MHz = 56.001 / 57.001 / 58.000 / 60.000 / 61.001 x 27.008 MHz; each line is weak and wobbles
  (appears as two lines ~36 kHz apart / short "bursts").
- **False positive "Iridium-like TDMA burst" (1/3) at 1620.48 MHz** = the 60th harmonic of that comb (no 8.28 ms
  bursts on a 90 ms frame, no 25 ksym/s QPSK; a < 15 kHz line wandering in frequency). -> Iridium classifier fix
  candidate (require burst length/frame timing and the 41.67 kHz channel raster, or reject comb lines).
- 1600.06 strong clock line (see 1588), 1640.0 integer line, 1644.45 weak.

### 1660 / 1696 MHz (1632–1724) — meteorological satellites 1675–1710
- No satellite downlinks. **1673–1687 MHz: impulsive broadband energy** (real per dither, 4–10 dB, kurtosis 6),
  |x|^2 lines at 96.0 / 192.0 kHz = 61.44 MHz / 640 -> capture-locked USB-link self-noise radiating at a fixed RF
  (same 96 kHz line as 1439.19). Engine: mostly below 6 dB; one 1.2 MHz `unknown` (1/3).
- 1701.52 MHz real, ~0.5 MHz, 5.6 dB, 3/3: FM-baseband 47.3 Hz wobble (as 1485.42) -> same EMI family; not MetOp
  AHRPT (1701.3 MHz, 2.33 Msym/s QPSK).
- 1714.9–1716.1 (660–720 kHz) and 1722.67 MHz (1.2 MHz, 16.9 dB), 1/3 each: B3 uplink (1710–1785) handset
  bursts. 1658.88 = 27 x 61.44 (labelled), 1640 / 1720 integer lines; 1647.5 = 61 x 27.008 comb.

### 1732 MHz (1704–1760) — B3 uplink (1710–1785) + H3 image of 5 GHz Wi-Fi
- ~60 bursty 0.4–3 MHz `unknown` signals (RMS 10–27 varying). Two sources:
  - **H3 image of an in-room 5 GHz Wi-Fi AP**, decoded through the image (`image_decode.py 1732 --n 3`, conj IQ at
    3 LO = 5196 MHz): **802.11 OFDM ch 36 (5180 MHz) beacon, BSSID c4:27:28:fa:d1:8f (success, 2/3); ch 40/44 VHT
    and NON_HT frames, client 74:0e:a4:9d:0e:49**. Image position F = 4 LO - source: ch 36 -> 1738–1758 MHz,
    ch 40 -> 1718–1738 MHz. (This AP is also why the ADC clips at 5172–5235 MHz.)
  - real B3 handset uplink bursts (1708–1718 MHz and elsewhere).
- **No signal carries `possible_lo_harmonic_image`**: the band-plan warning only knows cellular sources (STRONG_CLASSES)
  -> the image-decode confirmation (see 1156) would catch this; a 5 GHz Wi-Fi band-plan entry would over-flag.
- The dither map is unreliable for bursty signals (correlations 0.1–0.8 with mixed verdicts).

### 1768 MHz (1740–1796) — top of B3 uplink
- Few handset bursts (1742.3, 1750.7). The 27.008 MHz comb explains the wobbling lines: **1485.42 / 1701.52 /
  1755.52 / 1782.54 MHz = 55 / 63 / 65 / 66 x 27.008 MHz** (300–500 kHz wide, 47.3 Hz FM wobble, real, 3/3) ->
  one external oscillator (spread/wobbling 27 MHz clock), not emissions. 1781.76 = 29 x 61.44 (own).
  (1439.19 = 53.29 x 27.008: not on the comb; stays unidentified.)

### 1804 MHz (1776–1832) — B3 downlink starts at 1805 (strong: RMS ~550, clipping 1e-5)
- **LTE B3 FDD PCI 472, EARFCN 1251 (1810.1 MHz), 10 MHz (1805.1–1815.1)**: Engine partial (PCI) 2/3.
- **GSM-1800 ARFCN 566 (1816.0 MHz)**: FCCH 3/3 (purity 0.94–0.95, +7.7–7.8 ppm), SCH training 59/64 but CRC fails
  (MLSE gap). ARFCN 568 "FCCH" (1/3, purity 0.56, -9.8 ppm) = false (as at 976).
- **Engine miss: 1817.5–1832+ MHz, 40–42 dB above the floor, flat (a wide carrier running past the band edge), is not
  reported** — only 1816.7 / 1818.2 MHz fragments (8 dB). With ~27 of 56 MHz occupied the detector's floor
  estimate is likely pulled up. (Checked at LO 1840.)
- 1773.9 MHz, 1.3 MHz, 22.8 dB, 3/3: steady carrier in the B3 uplink band (seen at the edge; unidentified).

### 1840 MHz (1812–1868) — B3 downlink
- Occupancy: LTE 20 MHz 1817.5–1835.5, LTE 10 MHz 1836.5–1846, GSM 1815–1817, **1846–1868 empty**.
- **LTE B3 FDD PCI 110, EARFCN 1565 (1841.5 MHz), 10 MHz**: Engine partial 2/3; bench `lte_probe` decisive 3/3
  (SSS z 8.2–9.5 vs 4.6–5.2).
- **LTE B3 FDD PCI 101, EARFCN 1415 (1826.5 MHz), 20 MHz: bench decisive 3/3 (SSS z 6.2–6.3 vs 3.4–3.7, PSS z
  75–100); Engine: not detected at all** (only "weak detections" at ~0 dB across 1817–1835). Same PCI 101 as the
  B8 cell at 949.0 MHz (co-sited). PCI 472 (1810.1) seen at LO 1804.
- **Root cause (library, detect/noise.py)**: without a sweep NoiseModel (single block, e.g. `Engine.analyze_iq`) the
  floor is the 3 MHz morphological opening of the per-bin 20th-percentile profile, which by design "misses always-on
  signals wider than it": a flat 18 MHz carrier becomes the floor, so no candidate exists and the LTE receiver
  (which searches detector candidates) never tries it. Fix candidates: a single-block floor that uses the
  band-wide low percentile (here 1846–1868 is empty) / a stored front-end noise shape per device+gain (this live
  sweep's quiet LOs give it), and/or LTE/NR raster search inside wide flat regions.
- **Library bug (detect/artifacts.py clipping)**: for input normalised to [-1, 1] the threshold is
  0.995 * (fs_ - 1) = 0 -> clipped_fraction = 1.0 and ADC_CLIPPING on every normalised block (raw-count input,
  as signal-vault sends, is fine). Seen on the bench runs that pass x / 2048.

### 1876 MHz (1848–1904)
- Empty apart from the edges (1845.9 = top of PCI 110's carrier; 1905.7 truncated, see 1912): the upper B3
  downlink (1846–1880) is unused here. Dither "H7/H9" segments are noise-level (< 5 dB).

### 1912 MHz (1884–1940)
- No emissions. 1920.01 MHz, 840 kHz, 10 dB, 3/3 = 48 x 40 MHz (external 40 MHz comb, spread like 1600). The
  1905.7 MHz edge signal of LO 1876 is absent -> alias. No B1 uplink activity in these blocks.

### 1948 MHz (1920–1976) — B1 uplink
- No emissions (no handset uplink in these blocks). 1966.08 = 32 x 61.44 (own, labelled 3/3).

### 1984 MHz (1956–2012)
- No emissions. 1966.08 = 32 x 61.44 (own). B1 downlink starts at 2110 MHz (next session).

## Summary 896–2012 MHz (2026-10-10, LO 904 … 1984, 36 MHz steps, gain 40, --ppm 10)
Decoded/identified: LTE B8 PCI 96, 101; LTE B3 PCI 472, 110, 101 (bench only); GSM-900 BCCH ARFCN 97/100/102/106
(BSIC 47 on 106), GSM-1800 ARFCN 566; GSM handset uplink ARFCN 48 TSC 3; ADS-B 4 aircraft (1 s dwell); through
the LO 3rd-harmonic image: NR n78 PCI 972 (SFN 594/632/668) and 5 GHz Wi-Fi ch 36 BSSID c4:27:28:fa:d1:8f.
Everything else: own clocks (k x 61.44, 33.333 / DDR 1066.66), external combs (27.008 MHz x 55–66 with 47 Hz
wobble, 40 MHz x 40/48, 12 MHz family), USB self-noise (96 kHz-locked, 1673–1687 MHz), 3 unidentified weak lines.
Quiet: 984–1564 MHz apart from the n78 images (no DME, SSR, radar, GNSS, Inmarsat above the noise indoors).

## Fix list (library candidates from 896–2012 MHz; bench -> Python -> C++ with parity)
1. **Single-block floor misses wide always-on carriers** (detect/noise.py opening 3 MHz): LTE B3 PCI 101 (20 MHz,
   40 dB) not detected -> band-wide low-percentile floor / stored front-end shape; LTE/NR raster search in wide
   flat regions.
2. **GSM FCCH needs a CFO search** (fixed ppm off by 2.4 kHz hides all BCCH carriers) and an FCCH consistency gate
   (residual ppm within +-2 of the block, purity >= 0.8 or SCH training match) against false FCCH on neighbours;
   MLSE for SCH (SCH training 57–63/64 but CRC fails on most).
3. **GSM uplink burst detector** (GMSK 577 us on the 4.615 ms grid, 880–915 / 1710–1785 MHz, TSC -> BCC).
4. **False positives**: analog_vtx on a capture-locked fs/512 line at the LO (success at -18 dB SNR); Iridium-like
   on a comb line; DME-like on a narrow comb at 0.5 confidence.
5. **Image confirmation by decoding** (conj IQ at 3 LO / IQ at 5 LO for flagged or broadband energy) -> label the
   image with its decoded source; Wi-Fi images are not flagged today.
6. **Receiver-spur label depends on the detected width / floor** (k x 61.44 labelled in 2/3 snapshots).
7. **clipping() on normalised input** returns clipped_fraction 1.0 (threshold 0.995 x (1 - 1)).

## Correction (GSM FCCH, same day, after the fix list)
The "false FCCH" calls above (ARFCN 44–54 / 98 / 101 / 109 at LO 940 and 976, 568 at 1804) were wrong: their odd ppm
came from `find_fcch`'s lag-1 phase CFO, which residual adjacent-carrier power in the 220 kHz channel biases by
kHz. With the CFO taken from the interpolated FFT peak (`gsm_fcch_cfo.py`), every FCCH hit on LIVE_940MHz agrees on
+9.56..+9.80 ppm (13–16 carriers: ARFCN 11, 12, 44, 46, 47, 48, 50, 52, 54, 98, 100, 101, 102, 104, 106, 109) ->
they are real BCCH carriers (several co-sited: same FCCH instant). ARFCN 48 is a BCCH carrier, not only a TCH.
With the accurate CFO, ARFCN 102's SCH passes its CRC: BSIC x7 (NCC 4, BCC 7, as 106 at LO 940). Noise-only
channels peak at purity 0.32–0.39.

## 2012–3000 MHz (2026-10-15 ~23:10–23:50 local, LO 2020 … 2992, 36 MHz steps, gain 40 (30 at 2452), --ppm 10)
Night-time capture: capacity carriers may be off (B1 EARFCN 465 / 515, decoded via the 433 MHz image on 10-10, absent).

### 2020–2092 MHz (1992–2120) — space operations / B1 uplink top
- Quiet. 2033.335 = 61 x 33.333 (own). Weak narrow lines 2029.49 (1/3), 2071.89 (in LO 2056 and 2092: real, ~20 kHz),
  2088.89 (2/3): unidentified.

### 2128 / 2146 / 2164 MHz (2100–2192) — LTE B1 downlink
- **Decoded: PCI 103 EARFCN 315 (2141.5, 15 MHz, +40 dB), PCI 415 EARFCN 540 (2164.0, 10 MHz), PCI 426 EARFCN 565
  (2166.5)** - 415 and 426 overlap in frequency: different sites of one operator with different bandwidths (likely).
- 2131.6 MHz UMTS-like 4.4 MHz (+11 dB): WCDMA not decoded (no CPICH receiver).
- **2149–2158 MHz, +13–16 dB, not reported at all**: the energy region 2133.8–2158.5 is one detection; the LTE raster
  search finds only PCI 103. Bench PSS scan: weak cell at **2151.5 MHz (EARFCN 415), N_ID2 1, z 14–21 in 3/3
  snapshots**, SSS not decisive -> correctly not reported (PSS-only needs z >= 40). 2154–2158 unresolved.
- LO 2128 and 2164 see the 2133–2159 block truncated at the edge ("LTE-like ... truncated"): a dedicated LO was needed.

### 2200 / 2236 / 2272 MHz (2172–2300)
- Quiet. 2211.84 = 36 x 61.44 (own, labelled correctly). 2194.78 narrow 2/3. At LO 2272 the IQ image of the B40 edge
  (2244.7) is reported twice: once `iq_image`, once `unknown`.

### 2308 / 2344 / 2380 MHz (2280–2408) — LTE B40 TDD
- **Decoded: PCI 103 EARFCN 39150 (2350.0; same PCI as the B1 cell: same site/sector), PCI 350 EARFCN 39294 (2364.4,
  10 MHz).**
- 2309.8 MHz 20 MHz TDD carrier: PSS N_ID2 1 (z high), **SSS ambiguous** -> PCI not resolved; 2329.5 LTE-like TDD
  (truncated at LO 2308, not resolved at LO 2344).

### 2416 / 2452 / 2488 MHz (2388–2516) — 2.4 GHz ISM
- **Wi-Fi decoded:** 802.11b ch 2 "Airtel_upen_1325" (b4:a7:c6:45:21:d1), ch 6 "Airtel_Vini_1802", ch 10
  "Whitehouse"/"JioFiber-BaLA3_5G", ch 13 "Tata sangwan_2.4Ghz" (28:25:5f:99:f7:bb); OFDM ch 1 "Airtel_JSINGH37",
  ch 6 "Rekha"/"Advik" (6e:22:f7:f0:33:ac, 72:22:f7:f0:33:ac); ch 11 PHY headers only.
- **BLE advertising decoded (CRC) on 37/38/39** (2402, 2426, 2480).
- **Mislabel:** ~15 GFSK ~1 MHz bursts on data channels 2451–2478 MHz labelled "Bluetooth LE advertising-like":
  advertising exists only on 2402/2426/2480 -> BLE connection events or BR/EDR hops ("Bluetooth-like").
- Duplicates: `unknown` 1.3–3.8 MHz entries on top of decoded BLE adverts (2402.02, 2480.14 at +46 dB).
- No drone links (night).

### 2524 / 2560 / 2596 / 2632 MHz (2496–2660) — B41
- **Decoded: LTE B41 TDD PCI 218 EARFCN 40940 (2625.0).** At LO 2596 the dither map calls its truncated edge an
  "H3 image" (corr 0.62): the dither test is unreliable for signals cut by the band edge.

### 2668 MHz (2640–2696) — own clock comb
- 2666.68 MHz +30 dB with sidebands every 1.92 MHz (2657.08 … 2674.36): stays put under a 2 MHz LO dither ->
  RF, and **2666.67 = 80 x 33.333 = 5 x DDR 533.33 MHz: the Pluto's own clock**, modulated at 61.44/32 MHz. A second
  comb at LO + k x 1.92 MHz moves with the LO (baseband-locked self-noise). Library: `unknown` (the sidebands make
  the detection 1 MHz wide, above the 25 kHz receiver-spur limit).

### 2704–2992 MHz (2676–3020) — S-band radar
- Quiet except 2949.12 = 48 x 61.44 (own, labelled) and:
- **S-band ATC radar at 2888.57 MHz** (`dwell.py`, 3 x 2.1 s at 4 MSps; beam caught once, 78 ms): staggered PRI
  1117–1170 us; each group = ~1 us short pulse + **69 us LFM chirp (-450 -> +335 kHz, ~0.8 MHz)** = pulse
  compression (solid-state ASR, short pulse for near range). Rotation > 2 s (missed in 2 of 3 dwells).
  Library (LO 2884, 1/3 snapshots): `unknown`, 15 bursts, duty 0.19, band prior "S-band ATC / weather radar", but
  "no rule matched" and the chirp's autocorrelation read as "OFDM, SCS 140.6 kHz".

## Summary 2012–3000 MHz
Decoded: LTE B1 PCI 103 / 415 / 426; LTE B40 PCI 103 / 350; LTE B41 PCI 218; Wi-Fi SSIDs/BSSIDs on ch 1/2/6/10/13;
BLE adverts (CRC). Identified (bench): S-band ASR radar 2888.57 MHz (staggered PRI, 69 us LFM); weak B1 cell at
2151.5 MHz (PSS only). Not decoded: UMTS 2131.6, B40 2309.8 (SSS ambiguous) and 2329.5. Own lines: k x 61.44,
33.333 family incl. the 2666.67 DDR comb. No drones (night).

## Fix list (2012–3000 MHz)
Done 2026-10-16 (Python + C++, parity): 1 (radar: `features.pulse_train`, "Pulsed radar (S-band ATC)-like" with
PRI / LFM evidence, testdata `LIVE_SBAND_RADAR_2884MHz`) and 2 (labels: "Bluetooth (BLE/BR)-like" off the
advertising channels, testdata `LIVE_ISM_2452MHz`; the duplicate `unknown` entries remain).
1. **Radar classification:** pulse groups on a staggered ~1.1 ms PRI with an LFM long pulse -> "S-band ATC radar-like"
   (2700–2900 MHz prior); an LFM chirp must not read as OFDM; burst_period_s for staggered PRI.
2. **Bluetooth labels:** "advertising" only on 2402/2426/2480 MHz; 1 MHz GFSK hops elsewhere -> "Bluetooth-like"
   (BLE data / BR/EDR); absorb the duplicate `unknown` around decoded adverts.
3. **Own-clock lines with sidebands** (2666.67 DDR comb): a receiver-spur check by centre frequency, not width.
4. **LTE TDD SSS** (B40 2309.8 MHz ambiguous while 2350/2364 resolve): check SSS placement/averaging for TDD.
5. **Weak carrier next to a strong one in one energy region** (B1 2149–2158 MHz): per-carrier sub-regions
   (flat-profile split) so weak carriers get their own bandwidth/label even when the PSS is too weak.
6. Duplicate image reports (`iq_image` + `unknown` for the same mirror, LO 2272).

## 3000–4000 MHz (2026-10-16 ~00:30–00:55 local, LO 3020 … 3992, 36 MHz steps, gain 40, --ppm 10)
No clipping anywhere at gain 40 (n78 included, unlike the 2026-10-04 sweep).

### 3020 / 3056 MHz — pulsed emitter at 3036–3047 MHz
- Library: "Pulsed radar-like" 0.65 (1/3 snapshots at both LOs: same RF, so not an image), band prior maritime /
  radiolocation radar.
- Bench (`dwell.py` 4 x 0.84 s at 10 MSps + 0.14 s at 61.44 MSps; quiet control dwell at 3100 MHz: no pulses):
  pulses of **10–15 us, PRI exactly 4.000 ms (250 Hz)**, energy **3036–3047 MHz (~10 MHz, centre ~3041.5)**,
  present in every dwell over several minutes (peak 16–42 dB: no clean rotation pattern seen). Type not
  identified (radar-like by the pulse structure; no LFM resolved inside the pulses).
- Library flaw: the evidence still appends the autocorrelation reading "-> OFDM, SCS ..." for the pulsed signal.

### 3200 MHz — own DDR clock (6 x 533.33 MHz)
- 3200.01 MHz with ~80 kHz sidebands, fixed RF under LO 3190…3210: **6 x DDR 533.33 = 3200.0** (same family as the
  2666.67 comb). Library: `unknown` 2.5 MHz (the sidebands make it too wide for the receiver-spur check).

### 3272 MHz — LO-dependent own-clock beat (not RF)
- Two constant-envelope lines +28–30 dB at 3256.01 / 3288.02 MHz = **LO -/+ 3 x (LO - 3266.67)** (LO 3260…3274:
  pair at LO -/+ 20, 4, 10, 13, 16, 19, 22 MHz; absent at LO 3236 / 3276+). 3266.67 = 98 x 33.333 MHz: a beat of
  the LO with the Pluto's own clock harmonic, symmetric about the LO (real-valued), floor RMS 7 -> 16. Library:
  two `unknown` (and the dither map calls one an "H3 image of 9799.7 MHz").

### 3308 … 3596 MHz — 5G NR n78 (TDD)
- **Decoded: NR n78 PCI 102, SSB 3511.20 MHz (NR-ARFCN 634080), MIB CRC ok in 5/6 snapshots (LO 3488, 3524).**
- **PCI 531, SSB 3331.20 MHz (NR-ARFCN 622080)**: PSS/SSS only, PBCH CRC failed (1/3 at LO 3344; weak, 4 dB).
- Wide TDD carriers ~3380–3600 MHz: "5G NR-like TDD n78" without PCI, mostly truncated at the LO edges, and
  **fragmented**: one carrier becomes several `unknown` pieces (0.7–42 MHz) next to the NR-like entry.
- PCI 972 (SSB 3449.28, decoded through the H3 image at LO 1156 on 2026-10-10): not found directly at LO 3452
  ("5G NR-like" 54 MHz, truncated).

### 3632–3992 MHz — quiet
- Only 3932.16 = 64 x 61.44 MHz (own, labelled receiver spur correctly).

## Summary 3000–4000 MHz
Decoded: NR n78 PCI 102 (MIB). Partial: NR n78 PCI 531 (PSS/SSS). Identified (bench): pulsed emitter 3036–3047 MHz
(10–15 us, PRI 4.000 ms); own lines 3200.0 (DDR) and the LO-dependent beat near LO 3266.67. Nothing above 3600 MHz.

## Fix list (3000–4000 MHz)
7. **Own-clock lines, LO-dependent beat:** a symmetric pair at LO -/+ 3 x (LO - k x 33.333 MHz) -> receiver spur;
   3200.0 (6 x 533.33) with sidebands -> receiver spur (with item 3).
8. **Pulsed signals:** no autocorrelation "OFDM, SCS" reading in the evidence once a pulse train is measured.
9. **Wide NR TDD carriers fragmented** into `unknown` pieces: merge into the carrier (one carrier = one signal);
   SSB search over the whole carrier when the LO cuts it.

## Re-check of the known signals on fresh captures (2026-10-16 ~00:55–01:15 local, `--ppm 10`)
Found and decoded again: NR n28 PCI 75 (MIB); NR n78 PCI 591 at SSB 3449.28 MHz (MIB; PCI 972 was decoded on the same
NR-ARFCN through the image on 10-10) and PCI 102; LTE B1 PCI 415 (+ PSS-only 2141.5 / 2166.5), B3 PCI 350 / 230 / 76,
B5 PCI 6, B8 PCI 278 / 348, B40 PCI 103 / 350 / **196** (2329.8 MHz, unresolved on 10-15), B41 PCI 218; GSM-900 FCCH on
13 ARFCNs; DMR colour code 0 (858.05) and NXDN (856.82 / 861.57); ADS-B (80052B, 885281); Wi-Fi SSIDs and BSSIDs;
BLE adverts. Not seen: LTE B28 PCI 351 (weak), the S-band ATC radar (beam missed / off), TETRA MNC 4040 carriers.
- **TETRA MNC 2024 colour code 1**: a second network on 395.159 / 395.409 / 395.905 MHz, 81 / 61 / 33 CRC-valid BSCH
  with consistent slot / frame / multiframe counters (`testdata/live/LIVE_TETRA_395.28MHz`). MNC 1220 colour code 1 at
  390.079 MHz was decoded once (SNR 3.4 dB) and not again in two later dwells: unconfirmed.
- 860.821 MHz (DMR on 10-10) carried analog FM with CTCSS 171.3 Hz: a mixed-mode repeater.
- **False label: "sik like" on GSM-900 traffic bursts** (957.0 / 958.0 MHz, 554 us constant-envelope GMSK on the
  200 kHz raster, measured 154 kBd / 56 kHz deviation; no FCCH, so not claimed by `gsm`).
- Duplicate: two "BLE advertising (ADV_IND)" entries for one 2480 MHz advert (0.97 and 0.5).

## Fixes 2026-10-16 (Python + C++, unit tests, parity)
- Fix list 7 (part): none yet for own-clock lines (3200.0, 2666.67, the 3266.67 beat).
- **8 done:** an autocorrelation lag longer than a quarter of the median pulse width is not OFDM
  (`features.measure.OFDM_MIN_LAGS`).
- **9 (pulsed part) done:** narrower detections of one pulse train (PRI within 1 %, pulse starts within 5 us for half
  of the pulses, spans <= 2 MHz apart) and the burst-only slivers on its pulses -> one signal
  (`detect.merge.merge_pulse_trains`): the 3036–3047 MHz emitter is one "Pulsed radar-like" signal. The n78
  fragmentation (OFDM, not pulsed) remains.
- **GSM bursts no longer "sik like"**: `rc_link` drops 1–8 GSM time slots on the 200 kHz raster in GSM bands
  (`testdata/live/LIVE_GSM_958MHz_fsk`).
- Python regression (all 30 earlier inputs): no label changes; C++ parity identical on the new and affected inputs.
- Not changed: the TETRA decoder (no false positive shown).

## 4000–5000 MHz (2026-10-16 day, LO 4020 … 4992, 36 MHz steps, gain 40, --ppm 10, laptop kept awake)
Quiet below 4800 MHz: no radio altimeter (4.2–4.4 GHz, airborne only), no n79, no fixed links. Own lines labelled
correctly: 4066.0 / 4118.0 / 4868.0 / 4940.0 (reference-locked), 4177.92 (68 x 61.44), 4669.44 (76 x 61.44).

### 4871.7–4888.3 MHz — Mimosa fixed-wireless link (802.11 PHY), decoded on the bench
- Fixed RF under LO 4850…4900 (10 steps); 3.2 us CP autocorrelation, bursty envelope (p10 -1 / p90 +13 dB).
- `ofdm_rx.WifiDetector(scan_step_hz=1e6)` (blind 1 MHz grid) + `plugin_ofdm._records`: **238 NON_HT 6 Mb/s PPDUs
  at 4880.0 MHz (20 MHz), 230 FCS-valid** in 3 snapshots: 200 Block-Ack (56 B, TID 8) and 30 PS-Poll (20 B),
  **TA 20:b5:c7:02:d8:ca -> RA/BSSID 20:b5:c6:46:e8:3e (Mimosa Networks OUI)**, one frame every ~390–490 us: Mimosa's
  TDMA polling on an 802.11 PHY (point-to-point / multipoint backhaul; not a drone).
- **Library: missed.** `wifi_ofdm` runs only in the Wi-Fi bands and on the standard channel grid (4880 is neither), so
  the link is 2–3 `unknown` pieces of 4–5 MHz; `wifi_ofdm_scaled` rejects it (its bands are fragments: "occupied
  width 4.95 MHz does not fit 52-56 subcarriers").

### Own lines 4800–4960 MHz (LO steps 4840…4930)
- 4900.0 (+34 dB, 147 x 33.333) and 4800.03 (9 x 533.33 DDR): fixed RF, `unknown` (360 kHz wide with sidebands).
- 4915.2 (80 x 61.44, 300 kHz wide) and 4853.76 (79 x 61.44, at the LO 4884 band edge -> "truncated"): `unknown`.
- A comb at LO +-10 / +-20 MHz at some LOs (4890, 4910) and a floor rise (RMS 10 -> 15–20) at LO 4850–4920.
- 4423.68 (72 x 61.44) at the LO 4452 band edge: `unknown ... truncated` instead of a receiver spur.

## Summary 4000–5000 MHz
Decoded (bench): Mimosa 802.11-PHY link at 4880 MHz (MACs, 230 FCS-valid frames). Otherwise only own lines.

## Fix list (4000–5000 MHz)
10. **802.11 outside the Wi-Fi bands / off the channel grid:** run `wifi_ofdm` with the blind 1 MHz centre scan where a
    wide bursty signal shows a 3.2 us CP correlation (or add the 4.9 GHz band with the 802.11j grid 4000 + 5 x ch);
    expected: "Wi-Fi 802.11 OFDM ... Mimosa" at 4880 MHz.
11. Own lines (with item 7): sample-clock harmonics wider than 25 kHz (4915.2) and at a band edge (4853.76, 4423.68)
    -> receiver spur, checked before the "truncated" reason; 33.333 / 533.33 family with sidebands (4800, 4900).

## 5000–6000 MHz (2026-10-16 day, LO 5020 … 5992, gain 40, --ppm 10; paused at 5668 and resumed)
Quiet 4992–5490 MHz (no Wi-Fi on ch 36–64 at this antenna position, no weather radar). Own lines, `unknown`:
5100.0, 5200.0, 5400.0 (33.333 MHz family), 5320.0 / 5480.0 / 5640.0 (a 160 MHz-spaced comb); labelled spurs:
5406.72, 5652.48, 5898.24 (k x 61.44), 5864.0, 5904.0 (reference-locked).

### Wi-Fi (all decoded by the library)
- **5510.0 MHz, SSID "BI-Armory", BSSID 78:9a:18:bf:e8:4a** (40 MHz ch 102 centre: found by the off-grid scan).
- **ch 149 / 153 / 157: SSID "Armory" (+ a hidden SSID), BSSID 8c:86:dd:b7:4e:ab / ce:86:dd:b7:4e:ab**, 802.11ax
  (HE) + non-HT, stations 08:f9:7e:99:5e:dd, ac:f2:3c:ca:82:a1, 0a:39:90:e7:fb:70, bc:c7:46:e8:a8:4c.
- **Fragmentation:** one ~600 us wideband PPDU on ch 149 was reported as ~10 `unknown` pieces of 0.6–3 MHz at
  5739–5754 MHz (all with the same two bursts, overlapping the decoded PPDUs in time; LO 5755 snap00), and similar
  pieces at 5782–5814 MHz. **False label: "Pulsed radar-like" at 5495.65 MHz** = edge fragments of the 5510 MHz
  Wi-Fi (three 0.84 s dwells: irregular 1–15 us bursts, no PRI).

### 5945 MHz — analog FPV video (decoded)
- **NTSC analog video on band E ch 8 (5945.24 MHz), demodulated to a picture** (LO 5920, 3/3, buffer RMS 86); gone
  at LO 5956 a minute later (transmitter switched off / moved). Testdata already covers NTSC (Analog_20260804).

### Re-check captures added
`testdata/live/LIVE_NR_n78_3488MHz` (NR n78 PCI 445 MIB), `LIVE_LTE_B40_2344MHz` (B40 PCI 152 / 440, bench-confirmed),
`LIVE_MIMOSA_4880MHz` (off-grid 802.11).

## Fix list (5000–6000 MHz)
12. **One emitter = one signal for Wi-Fi:** energy pieces whose bursts coincide with decoded 802.11 PPDUs (and lie in
    their channel span) are absorbed into the Wi-Fi signal; this also removes the "Pulsed radar-like" fragment
    false positive in unii_2e (no radar label on a piece inside a decoded Wi-Fi span).
13. Own lines (with 7 / 11): the 160 MHz comb (5320 / 5480 / 5640) and 33.333 MHz multiples as receiver spurs.

## Fixes 2026-10-16 evening (Python + C++)
- **12 done:** `analyze._absorb_coincident` (>= 3 burst-only pieces in a receiver signal's span, same bursts +-20 us,
  >= 8 MHz together -> that emitter; LIVE_WIFI_5755MHz: 5 pieces into ch 149) and `rules.r_pulsed` needs a PRI in
  U-NII-2e (Wi-Fi fragments at 5495.65 MHz were "Pulsed radar-like").
- **10 done:** `wifi.ofdm_rx.offgrid_plan` (802.11 off the grid / outside the Wi-Fi bands; Mimosa 4880 / 5755 / 5790,
  BI-Armory 5510); trigger ~220 ms Python / ~70 ms C++ per block.
- Testdata added: LIVE_MIMOSA_4880MHz, LIVE_NR_n78_3488MHz (PCI 445 MIB), LIVE_LTE_B40_2344MHz (PCI 152 / 440),
  LIVE_WIFI_5755MHz.
- Open: own clock lines as receiver spurs (3 / 7 / 11 / 13), n78 carrier fragmentation (9), LTE TDD SSS (4), weak
  carrier beside a strong one (5), duplicates (6), pulse SNR for radars, latency (digital_vtx / nr_ssb / lte /
  wifi_ofdm_scaled take seconds per block in C++).

## 2440 MHz (2026-10-17 evening, LO 2440, gain 30, --ppm 10) — the 2.4 GHz "hopper"
- 94 µs 2 Mb/s GFSK packets (deviation ~±320 kHz as sliced) every 1.000 ms on exactly 2426.0 and 2450.0 MHz,
  20–40 ms per channel; two nodes 145 µs apart (sync 0x31A9D906 / 0x31A9D886: 3 address bits differ); a 3-bit field
  counting down mod 7 and two alternating payloads (idle traffic, repeats every 14 packets).
- CRC-16/CCITT (0x1021) found by a differential search over all 65,536 polynomials; it closes the first 124 bits
  after the sync word with an effective init 0xDD19 (86 / 88 packets). Plain ESB (init 0xFFFF over address + PCF +
  payload) and the XN297 scramble table do not fit: an ESB-family product of unknown make (wireless HID-like).
- Added to the library: `rc_link/esb.py` (+ C++), testdata `LIVE_ESB2M_2440MHz`. The other short bursts at
  2449–2465 MHz are Wi-Fi ch 10 pieces, those at 2420–2424 MHz their IQ images; 2432 MHz is a 1 Mb/s
  Bluetooth-like link; 2430.1 MHz is a second ESB-like device (sync 0xC7B6A60A).

## 890 MHz LoRa SF5 (2026-10-19) and ExpressLRS 2.4 GHz (2026-10-20)
- **890.03 MHz: LoRa SF5 / 500 kHz (SX126x)**, sync 0x12, 8-symbol preamble, ~0.2 ms carrier lead-in, one 6.67 ms
  packet every 112.0 ms, ~41 dB SNR. gr-lora_sdr has no SF5 frame (it decodes nothing). Frame measured here: SFD,
  2 extra symbols, first block CR 4/8 at the full SF (5/5 Hamming(8,4) codewords exact), then 15 blocks CR 4/5
  (75/75 parity-valid); **implicit header**: the first block is payload; 38 bytes + CRC-16, CRC-valid 12/12 in two
  dwells. Payload: bytes 11–14 = uint32 µs uptime (517 s, then 1.03 s after the user restarted the device), the rest
  near-constant (two message types). Likely the MicoAir telemetry radio (signal-vault `880MHz_MICOAIR_TELEMETRY`).
  Added: `lora_phy.decode_sf5_implicit` (+ C++), SF5 in the 500 kHz scan, testdata `LIVE_LORA_SF5_890MHz`, and the
  LoRa-only scan on 860–930 MHz (890 MHz is in B5 DL / B8 UL: before, the pipeline never ran rc_link there).
- **2440 MHz ExpressLRS "100Hz Full"** (SF7 / 812.5 kHz, inverted IQ, 13-byte packets every 10 ms, hopping over
  2412–2463 MHz on the ISM2G4 grid; one transmitter). The payload is not decodable: ELRS uses the SX1280
  long-interleaving coding rate `CR_LI_4_8`; the first block is reduced-rate (all symbols = 1 mod 4) but no
  standard/affine deinterleaver gives valid codewords. Only 2 distinct payloads in 11 snapshots. **Next:** a
  capture of a live link with moving sticks (dozens of distinct packets) -> find the interleaver from the parity
  checks (4-bit XOR sets that are zero in every packet), verify with the ELRS CRC-16 (poly 0x3D65, seed from UID).
  Added: testdata `LIVE_ELRS_2G4_2440MHz` (detection regression).
- Bench tool: `capture_testdata.py --burst-ms` waits for a buffer holding one whole burst (the pulse trigger fired on
  noise at 890 MHz).
- Regression (`bench/regress.py`, 43 inputs, old tree vs new): every signal identical except the intended
  "LoRa sf5 bw500k [success]" on `LIVE_LORA_SF5_890MHz`. The first version cost +42 s on the PMR dwell (SF5 false
  candidates on DMR/NXDN) -> SF5 stage-1 eta 40 + per-sub-band slice cache (see rc_link README).
