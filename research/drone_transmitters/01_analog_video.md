# Analog FPV video: standards, RF, decoders, identification cues (research 2026-10-21)

Status tags used below:
- **Read** sources: fetched and read.
- **[mem]**: from the researcher's memory of ITU-R BT.470 / BT.1700. The ITU PDFs could not be fetched.
- **[unverified]**: sources disagree, or the number comes only from vendor or forum listings.

Extracted source texts in this folder: `src_smpte170m.txt` (SMPTE 170M-2004, NTSC timing) and
`src_rtc6705_datasheet.txt` (Richwave RTC6705, the commonest VTX chip).

## 1. Composite standards after FM demodulation

| Parameter | NTSC-M | PAL-B/G/D/K/I | PAL-M (Brazil) | PAL-N (AR/PY/UY) |
|---|---|---|---|---|
| Line rate fH | 4.5 MHz/286 = **15 734.265 Hz** [S170] | **15 625 Hz** exactly | 15 734.265 [mem] | 15 625 [mem] |
| Field rate | 59.94 Hz | 50 Hz | 59.94 Hz | 50 Hz |
| Lines per frame / per field | 525 / 262.5 | 625 / 312.5 | 525 / 262.5 | 625 / 312.5 |
| Active lines | ~485 analog (480 digital) | 576 | 480 | 576 |
| Line period | 63.556 µs | 64 µs | 63.556 µs | 64 µs |
| H-sync width | **4.70 ±0.10 µs** | **4.7 ±0.2 µs** | 4.7 µs | 4.7 µs |
| Front porch | 1.5 ±0.1 µs | 1.65 (+0.4/−0.1) µs | ~1.5 µs | ~1.65 µs |
| Back porch / H-blank | sync edge to blanking end 9.20 µs; H-blank ~10.7 µs | back porch 5.7 ±0.2 µs; H-blank ~12 µs | ~10.9 µs | ~12 µs |
| Colour subcarrier fsc | 315/88 MHz = **3.579545 MHz** = 227.5·fH | **4.43361875 MHz** = 283.75·fH + 25 Hz | **3.57561149 MHz** = 227.25·fH | **3.58205625 MHz** = 229.25·fH + 25 Hz |
| Burst | **9 ±1 cycles**, starting 19 cycles (~5.3 µs) after the sync leading edge | **10 ±1 cycles** (2.25 µs), ~5.6 µs after the leading edge | 9 ±1 [mem] | 9 ±1 [mem] |
| Burst phase | fixed, 180° | **±45° line-by-line swing** (135°/225°) | PAL swing | PAL swing (+25 Hz offset) |
| V-sync | 6 pre-equalising + 6 broad + 6 post-equalising pulses at half-line spacing; equalising pulse 2.30 µs; serration 4.70 µs; field blanking 20 lines; no burst on the 9-line V interval | 5 + 5 + 5 pulses; equalising ~2.35 µs; broad ~27.3 µs; field blanking ~25 lines; Bruch blanking [mem] | NTSC-like 6/6/6 | PAL-like 5/5/5 |

**SECAM:** FM subcarriers Db 4.250 / Dr 4.40625 MHz, no burst. No FPV camera outputs SECAM.

**What separates the standards:**
- **525-line from 625-line** (NTSC/PAL-M vs PAL/PAL-N): line rate (0.70 % apart; measurable to ~1 ppm), field rate,
  lines per field, and the V-sync pulse count. **All of these work without colour.**
- **NTSC from PAL-M, and PAL from PAL-N:** chroma only.
  - fsc/fH: 227.5 vs 227.25 (3.93 kHz apart), and 283.75 vs 229.25.
  - PAL-family bursts swing ±45°; NTSC bursts don't.
- **Monochrome:** only "525/59.94" or "625/50" can be stated (NTSC-like / PAL-like). India uses PAL-B/G for
  broadcast, and FPV cameras are NTSC/PAL switchable.

## 2. VTX RF

**Modulation:** composite video FM-modulated onto the carrier. RTC6705 video input is 1 Vpp into 75 Ω. The datasheet
gives no video deviation and no video pre-emphasis; audio pre-emphasis only, corner ~12 kHz [unverified].

**Deviation (no manufacturer specification):**
- BATC gives a PAL FM-ATV channel width of ~18 MHz.
- By Carson's rule, with fm ~6.5 MHz, Δf ≈ ±3–4 MHz peak.
- The repo assumes a 4–10 MHz sync-to-white swing, which agrees. **Measure it rather than assume.**
- Polarity (sync tip at the low or high frequency) is not standardised: measure it per unit.

**Audio subcarriers:**
- RTC6705: two FM subcarriers at **6.0 and 6.5 MHz** (stereo, or mono on 6.5), **−25 dBc**, **±25 kHz** deviation.
- Most 5.8 GHz VTXs that carry audio use 6.5 MHz. 1.2/1.3 GHz VTXs: 6.0 / 6.5 MHz. 5.5 MHz also occurs.
- Many micro VTXs carry no audio.

**Synthesiser:** RTC6705 uses an 8 MHz crystal and a 20 kHz reference step, so channel centres are exact apart
from crystal ppm.

**5.8 GHz channels (MHz)** [OL, BF]:

| Band | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| A (Boscam A / TBS) | 5865 | 5845 | 5825 | 5805 | 5785 | 5765 | 5745 | 5725 |
| B | 5733 | 5752 | 5771 | 5790 | 5809 | 5828 | 5847 | 5866 |
| E | 5705 | 5685 | 5665 | 5645 | 5885 | 5905 | 5925 | 5945 |
| F (Fatshark / ImmersionRC Airwave) | 5740 | 5760 | 5780 | 5800 | 5820 | 5840 | 5860 | 5880 |
| R (Raceband) | 5658 | 5695 | 5732 | 5769 | 5806 | 5843 | 5880 | 5917 |
| D (Diatone / LowRace = Raceband − 296; Lumenier calls it "F"; **the repo's "L"**) | 5362 | 5399 | 5436 | 5473 | 5510 | 5547 | 5584 | 5621 |
| L (Lowband, 40 MHz step) | 5333 | 5373 | 5413 | 5453 | 5493 | 5533 | 5573 | 5613 |
| U | 5325 | 5348 | 5366 | 5384 | 5402 | 5420 | 5438 | 5456 |
| O | 5474 | 5492 | 5510 | 5528 | 5546 | 5564 | 5582 | 5600 |
| H | 5653 | 5693 | 5733 | 5773 | 5813 | 5853 | 5893 | 5933 |
| X [unverified] | 4990 | 5020 | 5050 | 5080 | 5110 | 5140 | 5170 | 5200 |

Band letters vary by vendor, so report all matching (band, channel) candidates, e.g. R7 = F8 = 5880.

**Other bands:**
- **1.2/1.3 GHz** [unverified]: 1080, 1120, 1160, 1200, 1240, 1280, 1320, 1360, plus 1258 MHz; hardware range 1060–1380.
- **2.4 GHz:** 2414 / 2432 / 2450 / 2468; 8-channel sets add 2370, 2390, 2490, 2510.
- **3.3 GHz** (3.1–3.5 GHz): no authoritative table; detect blind.

## 3. Decoders (all licences need checking; none permissive enough for C++)

| Tool | Licence | Use |
|---|---|---|
| PySDR ch. 25 FPV | CC BY-NC-SA | reference |
| SDRangel demodatv | GPL-3 | oracle |
| hacktv | GPL-3 | generator; transmit only into a dummy load |
| ld-decode | GPL-3 | oracle |
| SigDigger | GPL-3 | oracle |
| uAVD | closed | — |

## 4. Identification cues from IQ

**RF:**
- 15–20 MHz occupied bandwidth.
- A sync-tip dwell spike (sync is ≥ 7.4 % of each line).
- A line comb at fH in the spectrum.
- A channel-table match within crystal ppm.
- Audio sidebands at ±6.0 / 6.5 MHz, ~25 dB down.

**After FM demodulation:**
- **Line and field rates:** lines at k·fH with field-rate sidelines (59.94 / 50 Hz).
- **Sync and field structure:** a fold on the line period (4.7 µs sync); lines between V-sync (262.5 / 312.5);
  equalising-pulse count (6 / 5).
- **Chroma:**
  - NTSC interleaves at (n+½)·fH; PAL at (n±¼)·fH + 25 Hz.
  - Burst window ~5.3 µs (NTSC) / 5.6 µs (PAL) after the sync leading edge: constant phase (NTSC) vs ±45° swing (PAL).
  - PAL-M / PAL-N: fsc ≈ 3.576 / 3.582 MHz plus the swing.
- **Audio:** an FM subcarrier at 5.5 / 6.0 / 6.5 MHz, about 1.07 MHz above PAL chroma.

**Pitfalls:**
- Cameras switch between NTSC and PAL, so label the signal, not the drone.
- The 4:3 / 16:9 setting does not change timing.
- **OSD:** a MAX7456 OSD with no camera makes its own sync and a black field, which looks monochrome.
- "240p/288p" progressive sources have no half line.
- Farkas claims ±1 % line-timing variation (doubtful for crystal cameras): always confirm with field rate.
- **False positives at the line rate:**
  - the Pluto's capture-locked 15 625 Hz noise (8 MSps/512);
  - LTE (1 ms = 15.625 lines);
  - Wi-Fi (4 µs).

  Require a real V-sync and sync-tip dwell.
- DJI O3/O4, Walksnail and HDZero transmit only digital, with no analog fallback.
- **Images:** check Pluto 3rd/5th LO-harmonic images for off-table frequencies.

## Sources
- SMPTE ST 170M-2004: https://pub.smpte.org/pub/st170/st0170-2004_stable2010.pdf
- Wikipedia, PAL: https://en.wikipedia.org/wiki/PAL
- PySDR, FPV video chapter: https://pysdr.org/content/fpv_video.html
- RTC6705 datasheet: https://der-frickler.net/_media/modellbau/fpv/sky-rf/rtc6705-dst-001.pdf
- Oscar Liang, FPV channels: https://oscarliang.com/fpv-channels/
- Betaflight VTX docs: https://betaflight.com/docs/wiki/guides/current/VTX
- aldeid 5.8 GHz table: https://aldeid.com/wiki/Category:Drones/FPV
- BATC 5.6 GHz FPV ATV: https://wiki.batc.org.uk/images/9/9c/20190629_5.6_GHz_FPV_ATV.pdf
- FlyingTech FX796T VTX: https://www.flyingtech.co.uk/product/fxt-fx795t-2-5-8ghz-40ch-fpv-vtx-25mw-to-200-600mw/
- 1.3 GHz listings:
  - https://rotorvillage.ca/flywoo-1-3ghz-ultra-high-output-9-channel-vtx
  - https://rcdrone.top/es/products/1-2g-1-3ghz-1-6w-9ch-vtx-es
- 2.4 GHz listings:
  - https://www.flitetest.com/articles/Video_TX_RX_Frequency_Chart
  - https://blog.fpv-community.de/threads/verkaufe-fpv-groundstation-sender-2-4ghz.7632
- Farkas, Hadmérnök 2026: https://real.mtak.hu/248007/1/08_farkas_119-135_HM2026_2.pdf
- arXiv CNN FPV detection: https://arxiv.org/pdf/2607.16455
- PAL-M / PAL-N subcarrier patent: https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/5555305
- hacktv: https://dyn.manpages.debian.org/trixie/hacktv/hacktv.1.en.html
