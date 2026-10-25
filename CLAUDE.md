# spectrum_intelligence — notes for Claude

This is the user's primary development package. It detects, localises, identifies and decodes every signal in
IQ blocks (ADALM-Pluto and other SDRs, 300 MHz–6 GHz, India). `README.md` describes it; read it first.

Two implementations of the same library:
- `python/` (`spectrum_intel`): the development reference, used by the signal-vault GUI;
- `c++/`: the production library for edge devices (`c++/README.md`, `c++/PORTING.md`).
They share `schema/` (response JSON Schema, technology taxonomy, receiver requirements) and `testdata/`.

## Rules
- **Development order: bench → Python → C++.** A new signal is first detected and decoded on the bench. Only when
  it is verified there does it go into `python/` (validated on `testdata/`). Then the same tested, verified
  capability is ported to `c++/`. Code in `python/` and `c++/` is assumed tested and working: never add
  unverified work to either, and never let C++ carry something Python does not.
- **C++ follows Python.** Every verified Python change is ported to `c++/` and passes the parity check
  (`c++/tools/parity.py`: same signals, labels, decode status and decoded fields on testdata). Track each
  receiver's C++ status in `c++/README.md`.
- **Licences:** the C++ library must stay free of GPL/AGPL code. pyModeS (GPL-3) and py3gpp (AGPL-3) are test
  oracles only; implement from the standards. Third-party C++ code: pocketfft (BSD-3), nlohmann/json and doctest
  (MIT).
- **Response contract:** `schema/response.schema.json` is versioned (additive change = minor, rename/removal =
  major). Python (`spectrum_intel.api`) and C++ (`si::Engine`) must emit the same response. Expanding it when a
  new signal or decoded field needs it is allowed (additive, both implementations together).
- **Keep README.md current.** After every change to the library, update `README.md`: the signal table (what is
  detected and decoded, and how), the results and the limitations. Also update the receiver's
  `spectrum_intel/receivers/<name>/README.md`.
- **Stay modular and loosely coupled.**
  - A receiver's DSP imports only numpy/scipy and `spectrum_intel.core` (C++: only `core/`, `io/`).
  - `plugin.py` is its only adapter to the pipeline.
  - Nothing outside `decode/registry.py` imports receivers.
  - `python/tests/test_modularity.py` and `c++/tools/check_modularity.py` enforce this; keep them passing.
- **The machine has 8 GB RAM.** Run heavy jobs with `--workers 2` or fewer, under `nice -n 19`. Six workers froze
  the system before.
- **Never guess silently.** A label below 0.5 confidence is `unknown`, with evidence. An inferred label gets the
  "-like" suffix. Every report signal carries its `detection` and `decoded.steps` provenance.
- **Latency: as low as possible without losing accuracy.** The goal is every detection and decode in
  milliseconds per block. Every change is measured for run time (regression timings, C++ timings) as well as
  results. Reuse what the pipeline already computed (the block FFT / `Channelizer`, cached scans), gate expensive
  searches behind cheap triggers, and never trade a correct detection or decode for speed.
- **Add every new decodable signal.** A signal found on the bench that we can detect and decode but the library
  lacks (or has no regression capture for) is added: testdata in `testdata/live`, Python, C++ with parity, docs.
- **Every new signal gets testdata samples**, decodable or not: one block holding the whole burst/packet in
  `testdata/live/<NAME>` (SigMF ci16; `capture_testdata.py --trigger <MHz> --burst-ms <ms>`, or a crop of a dwell /
  signal-vault capture), plus a `testdata/README.md` row with the ground truth and the expected output, for later
  validation.
- **Prefer existing open-source decoders** (in C++ only permissively licensed ones; GPL/AGPL code only as a test
  oracle). Validate a decoder on real captures with known contents
  (`testdata/README.md`) and on synthetic signals before trusting its output on the sweep.

## Commands
```
.venv/bin/spectrum-intel run <path> --out out/full --workers 2 [--resume] [--freq-range LO:HI] [--disable R]
.venv/bin/spectrum-intel receivers
cd python && ../.venv/bin/python -m pytest -q tests
cd python/tests/receivers && ../../../.venv/bin/python test_det_<name>.py
cmake -S c++ -B c++/build -G Ninja -DCMAKE_BUILD_TYPE=Release && nice -n 19 ninja -C c++/build -j2
c++/build/si_tests && .venv/bin/python c++/tools/parity.py --cpp c++/build/si_analyze [--receivers R] [--inputs G]
```
`.venv` is radioconda Python 3.10 (`--system-site-packages`) with py3gpp, pyModeS, jsonschema, pytest and pyflakes;
`python/` is installed editable in compat mode (`pip install -e python --config-settings editable_mode=compat`), so
`python/` is on sys.path. `SPECTRUM_INTEL_THREADS` sets the threads per block (default min(4, cores)). Disk is
nearly full (~3 GB free): Release C++ builds only, delete scratch build directories.

## signal-vault (the viewer)
- **Where:** the acquisition app is a separate git repo at `../signal-vault`; this work is on branch
  `spectrum-intel-viewer`. Its **Signal Intelligence** tab shows this library's reports live, on replays and on
  files, and plays analog video.
- **How it connects:** the bridge is `signal_vault/intel/` (worker process + video thread) and the tab is
  `signal_vault/ui/intel.py`. They use only `analyze_iq` / `analyze`, `report.writer.build`,
  `decode.registry.receiver_names` / `plugins`, `receivers.analog_vtx.FieldDemodulator` and
  `receivers.voice_fm.VoiceDemodulator`. Keep those APIs and the report schema backward compatible.
- **Always the latest code:** the GUI imports the library from `python/` (editable install, or
  `../spectrum_intelligence/python` via `signal_vault/intel/source.py`). `source.Tracker` fingerprints the source
  before each analysis; after an edit the engine starts a fresh worker and re-imports the demodulators.
- **Environment:** signal-vault is installed editable into this repo's `.venv` (plus sigmf and
  opencv-python-headless). Run it with `.venv/bin/sv gui` or `.venv/bin/sv gui --replay <capture|folder>`.
- **Its tests:** `cd ../signal-vault && ../spectrum_intelligence/.venv/bin/python -m pytest -q tests`.
  `tests/test_intel.py` drives the tab end to end on `testdata/`.

## Data facts that matter
- **Production input:** IQ + sample rate + centre frequency + analog bandwidth (+ optional start time, gain,
  receiver ppm). Minimum 65,536 samples per analysis; recommended 2^21 (34 ms at 61.44 MSps). Smaller SDR blocks go
  through `StreamAnalyzer`. Without `receiver_ppm` a single block loses weak GSM/NR cells (testdata/cellular).
- **Sweep 2026-10-04:**
  - location: `/Users/srikanth/Documents/modules/signal-vault-data/captures/2026-10-04/pluto/SWEEP_300-6000MHz_SWEEP_Height0m_Distance0m_ControllerDist0m_Gain40dB_BW61.44MHz`;
  - 713 snapshots, 181 steps of 31.44 MHz, 4 repeats about 0.33 s apart, cf32 raw ADC counts, 1,998,848 samples;
  - full results are in `out/full/<sweep>/`.
- **Input is trusted.** The library assumes the data it receives is correct: the buffers were flushed and the LO
  matches the metadata at capture time. Do not add code that repairs input problems, such as stale-LO detection;
  they are fixed on the acquisition side (signal-vault). The 2026-10-04 sweep itself has stale buffers in 3 of the
  4 repeats per step and needs recapturing; its `out/full` results date from when the pipeline still corrected
  this.
- **Pluto TCXO:** about +7.4 ppm, measured from the LTE PSS at 940.1 MHz, and it drifts 1–2 ppm across the sweep
  (+9.5–10 ppm on 2026-10-10 afternoon: use `--ppm 10` / the LTE-measured value on the bench).
  `sweep.calibrate_ppm` measures it; receivers search a wider CFO grid.
- **Pluto front end (bench 2026-10-10, `bench/live_lo/FINDINGS.md`):** no preselector, so strong 1.8–2.6 GHz
  downlinks appear at 330–830 MHz through the mixer's 3rd (inverted) and 5th LO harmonics; a +1 MHz LO step moves
  them +4 / −4 MHz (`bench/live_lo/harmonic_map.py`). The USB link adds capture-locked self-noise (floor +4–5 dB
  ~12 ms into each 2^21 buffer, 8 kHz microframe modulation). Own clock lines: k × 61.44 MHz, 33.333 MHz (Zynq,
  DDR 533.33), integer-MHz reference lines. A retune leaves 1 stale buffer: destroy the rx buffer, discard 4.
- **Live bench tools:** `bench/live_lo/lo_step.py` (one LO: flushed capture + Engine + dither map),
  `capture_testdata.py` (SigMF into `testdata/live`), `dwell.py`, `triage.py`; `bench/regress.py run|diff`
  (Python reports on all testdata, before/after a change).
- **Clipping:** the ADC clips in strong cellular bands (770–960, 1808–1871, 2122–2185, 2280–2405, 3411–3537 and
  5172–5235 MHz).
- **Wrong labels:** some `drones_data` captures have wrong frequency labels (Phantom 4 +20 MHz, Herelink −10 MHz).
  The BLE label check (`receivers/rc_link/label_check.py`) catches this.
- **Origin:** the receivers came from the `wifi-decoding` research tree (`../wifi-decoding/*.py`), which is no
  longer used. Develop here.
- **DroneID:** the payload decoder is the external repo `../droneid-decoder-reference/src` (or
  `$SPECTRUM_INTEL_DRONEID_SRC`).

## Next work the user discussed
- continue the live sweep from 3000 MHz (n78, 5 GHz) with the bench tools (`step.sh`; 2012–3000 MHz done
  2026-10-15 at night, findings + fix list in `bench/live_lo/FINDINGS.md`: radar classification, Bluetooth labels,
  own-clock comb with sidebands, LTE TDD SSS, weak carrier beside a strong one);
  open items of the 896–2012 MHz fix list (`bench/live_lo/FINDINGS.md`): GSM uplink burst detector, image
  confirmation by decoding at 3 x LO / 5 x LO, receiver-spur label stable when the floor rises (done: single-block
  floor bridging, GSM FCCH CFO search + consensus, analog/Iridium/DME false positives, normalised clipping);
- TETRA SYSINFO (BNCH), DMR link control / vendor CSBKs, NXDN LICH;
- LTE PBCH/MIB and SIB1 (to identify the operator);
- a GSM MLSE equaliser for SCH/BSIC, then BCCH system information;
- NR SIB1;
- UMTS CPICH;
- DVB-T2 L1;
- grouping unknowns across snapshots, then an ML classifier for what remains unknown;
- acquisition fixes in signal-vault (buffer flush, per-band gain, longer dwell).
