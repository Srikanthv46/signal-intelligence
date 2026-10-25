# DJI O4 emitter identity: bench findings (2026-10-16)

## Question
One DJI drone hops between 2.4 and 5.8 GHz, switches 10/20(/40) MHz, and its received power swings with range.
Is it one drone? Can two drones of the same model be told apart?

## Data
All analysed at 61.44 MSps; the centre frequency comes from the filename; sidecar metadata ignored.
- `signal-vault-data/recordings/hopping`: one O4 drone at 06:58 UTC. Five simultaneous channels (2428, 2484, 5728,
  5784, 5840 MHz), 4 back-to-back 2,000,000-sample chunks. One 20 MHz link: 2415.5 MHz in chunks 0–2, then
  5823.5 MHz in chunk 3.
- `signal-vault-data/recordings/debug_detection`: the same drone at 07:10 UTC, in 10 MHz mode at 5730.5 MHz.
  It is weak, about 2–5 dB above the noise. Only chunk 0 is usable; the chunk-1 files are 1.8 ms long.
- Comparison drones:
  - `testdata/DJI_Mini_Pro_4_*`: O4, 10 MHz, 5770.5 MHz;
  - `tmpdata/D_20260609-151042` ch0: O4, 20 MHz, 5770.5 MHz.

## O4 structure (measured)
- **Slot structure.** A 0.5 ms slot has 13 symbols of 2192 samples plus one long-CP symbol of 2224 samples (CP 176),
  like NR at 30 kHz. A grid with 14 equal symbols of 2194.29 samples misplaces FFT windows by up to 16 samples.
- **Reference symbols.**
  - They sit at the odd slot positions 1–13, with pilots on alternate subcarriers (comb-2).
  - The set recurs at slot lags of 2 and of 5 (the 2.5 ms frame); nothing matches at a 1-slot lag.
  - Within a slot, the symbol at position r recurs 4 and 8 symbols later.
- **10 MHz mode** carries a 300-subcarrier slice of the 20 MHz sequence, offset by 150 subcarriers.
- **The sequence is the same on every O4 drone seen.**
  - Today's drone matches the Mini 4 Pro capture at p ≈ 1e-81.
  - It matches D_20260609 at p ≈ 3e-19.
  - So the sequence identifies "DJI O4" at any power, band or bandwidth, but **it does not identify the unit**.
  - Caveat: unless the Mini 4 Pro capture is the same drone (still to confirm).
- **Channel raster.** Channel centres lie on 0.5 MHz steps, all at x.5 MHz in these recordings (2415.5, 5823.5,
  5730.5, 5770.5).

## Method (`o4.py`)
Each link is measured over one contiguous stretch.
1. **Symbol grid:** the CP metric folded on the 0.5 ms slot.
2. **Fine CFO:** the CP phase, ±15 kHz.
3. **Sequence comparison:** comb-2 differential pilot spectra X(k+2)·X*(k), at unit modulus. They are unaffected by
   power, channel and timing, and are matched against another link over integer subcarrier shifts. The result is
   significant only when many symbols agree on one shift (binomial p).
4. **Symbol-clock ppm:** the timing drift of recurring reference symbols. It comes from a least-squares phase slope
   across the pilots, then a least-squares fit against time separation over every recurring pair. It needs no
   channel raster.
5. **Carrier:** f0 + CFO + 30 kHz × n.
   - n is chosen so that carrier / (1 + clock ppm) lands on the 0.5 MHz raster.
   - Sequence alignment against a reference link gives n independently: n_ref − shift for 20 MHz,
     n_ref − 150 − shift for 10 MHz.
   - Carrier ppm = (carrier − nominal) / nominal. **This carrier ppm is the crystal fingerprint.**

## Synthetic validation (`test_synth.py`, 8/8 pass)
| Check | Result |
|---|---|
| 20 MHz, 20 dB | clock +1.502 ± 0.002, carrier +1.4996 ppm (true +1.5) |
| 20 MHz, 10 dB | clock +1.507 ± 0.006, carrier +1.4989 ppm |
| 20 MHz, 5 dB | clock +1.529 ± 0.022, carrier +1.4983 ppm |
| 20 dB power ramp across the stretch | carrier −2.2995 ppm (true −2.3) |
| 10 MHz, 10 dB | carrier +1.4955 ppm |
| 10 MHz vs 20 MHz, same sequence | shift −150, p ≈ 1e-204 |
| Different sequence | no match |
| Two transmitters 0.2 ppm apart in one capture | measured difference +0.1988 ppm |

## Real recordings (`identity.py links_2026-10-16.json`)
| Link | Carrier ppm | Clock ppm | n (clock) | n (sequence) |
|---|---|---|---|---|
| Today 2.4 GHz, 20 MHz, chunk 0 | +1.459 | +2.39 ± 0.10 | 0 | 0 |
| Today 2.4 GHz, 20 MHz, chunk 1 | +1.462 | +1.79 ± 0.13 | 0 | 0 |
| Today 2.4 GHz, 20 MHz, chunk 2 | +1.423 | +2.68 ± 0.50 | 0 | 0 |
| Today 5.8 GHz, 20 MHz, after the hop | +1.536 | +0.79 ± 0.14 | −1 | −1 |
| Today 5.73 GHz, 10 MHz, 12 min later | +1.512 | too weak | — | +1 |
| Mini 4 Pro, other session and receiver | +12.42 | +11.98 ± 0.57 | +3 | +3 |

- **One crystal.** The carrier ppm of today's drone agrees within ±0.06 ppm across the band hop, the bandwidth change
  and 12 minutes. All numbers are relative to the recorder's clock.
- **Precision.** Per block it is about 0.01–0.03 ppm. Unit-to-unit spread in TCXOs is about 2–3 ppm.
- **Independent check of n and the raster.** Clock and sequence alignment agree on n for every link. The real-data
  clock scatters by about ±1 ppm, against ±0.002 on synthetic data: something not yet modelled (motion,
  multipath, recorder timing). It is still enough to fix n, since one subcarrier is 5–12 ppm, but it is not a
  fingerprint by itself.
- **The 2.4 vs 5.8 GHz difference** of about 0.08–0.1 ppm is larger than the per-block noise. Candidates: Doppler,
  a recorder channel offset, or CFO bias from Wi-Fi overlapping the 2.4 GHz link.

## What follows for identity
1. **Two same-model drones flying at once.** Each one transmits video more than half of every 2.5 ms frame. With all
   five channels recorded at the same time, two O4 links on air together means two drones. Their crystal ppm and
   their separate frame grids also separate them within a block.
2. **Same drone over time, or across hops and bandwidth changes.** Track the drone on carrier ppm (gate about ±0.2 ppm).
   A random second unit falls inside that gate about 8% of the time. Frame-timing continuity tightens this, but
   needs gap-free recording. Bandwidth, band and power are modes or nuisances, not identity.

## Open items
- **Recorder timing jumps.** The signal-vault recordings show symbol-grid jumps about every 1M samples and at every
  chunk boundary. Fit only between jumps; fix it in acquisition.
- **Real-data clock scatter.** Explain the ±1 ppm (vs ±0.002 synthetic).
- **Is the Mini 4 Pro capture the same drone as today's?** Determines whether "sequence common to all O4 drones" is
  proven.
- **Two-drone recording.** Fly two drones of one model: alone and together, cold start, forced bands, 10/20/40 MHz,
  near and far. This measures the real unit-to-unit spread and thermal drift, and gives a 40 MHz capture
  (`dji/detector.py` has no 40 MHz profile).
- **Then move to `python/`** (the `dji` receiver plus a cross-block tracker), followed by the C++ port and parity check.

## Short chunks of extracted packets (`chunk.py`, `chunk_study.py`)
**Setup.**
- Packets are cut with signal-vault's own `ddc_` at fs = 61.44e6, from the library's O4 bursts and boxes. For chunks
  1–2 the box is 3 MHz off-centre and too wide, like a real detector error.
- Two packet rates: box width (`PacketRate.BANDWIDTH`) and integer decimation to an O4-native rate
  (15.36 / 30.72 / 61.44 MSps).
- Chunk measurement:
  1. resample the chunk only, to the native rate;
  2. find symbols one by one from CP peaks;
  3. take the fine CFO from those CPs;
  4. find the integer bin with a template from another capture (never itself);
  5. carrier = mixing frequency + CFO + 30 kHz × n.

**Acquisition: the integer bin found from the chunk alone.** Totals for the four 20 MHz links (about 24 dB):

| Chunk | Samples at the box rate | Integer found | Wrong integer | ppm SD | 95% within |
|---|---|---|---|---|---|
| 3 symbols (110 µs) | 2,000–2,700 | 42–63% | **0** | 0.07–0.14 | 0.15–0.30 |
| 1 slot (0.5 ms) | 9,100–12,200 | 78–100% | **0** | 0.04–0.07 | 0.08–0.18 |
| 2 slots (1 ms) | 18,200–24,500 | 82–100% | **0** | 0.03–0.07 | 0.06–0.15 |

**Tracking: the integer bin known from the drone's track.** A drone's carrier is fixed on a channel and its ppm is
known to about ±0.1, so every chunk with symbols gives a ppm value.

| | 3 symbols | 1 slot | 2 slots |
|---|---|---|---|
| 20 MHz, about 24 dB: SD (95% within) | 0.07–0.11 (0.16–0.27) | 0.04–0.07 (0.07–0.17) | 0.03–0.07 (0.06–0.15) |
| 10 MHz, about 2–5 dB: SD (95% within) | 0.6 (1.4): unusable | 0.12–0.17 (0.27–0.31) | 0.11 (0.22) |

- **No wrong integers.** Across 1,036 chunks a failure is always a rejection (weak template match), never a wrong
  integer, so the method fails safe.
- **Box rate vs integer decimation:** identical results. Resampling only the chunk handles arbitrary box rates.
- **Design that follows.**
  - Acquire a new link with a 1-slot chunk; if it is rejected, take the next packet.
  - After that, short chunks (about 3 symbols at good SNR) are enough.
  - Average chunks over packets: N chunks give SD/√N.
  - Weak links (5 dB or less) need 1 slot or more, and several packets.
- **2048 samples at the box rate** is about 3 symbols at 20 MHz (18 MSps): fine for tracking. It is about 6 symbols at
  10 MHz (9 MSps). At 40 MHz (about 37 MSps) it is only 1.5 symbols, so use about 4,100 there; the 40 MHz numerology
  is still unverified.
- **Biases:** medians are ±0.04 ppm of the full-stretch value on the 20 MHz links. The weak 10 MHz link reads
  −0.12 to −0.25 ppm, partly because its full-stretch reference is itself noisy.
