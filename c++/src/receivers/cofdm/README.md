# cofdm (C++) — DVB-T-framed COFDM links (research capability)

Port of python `receivers/cofdm` (research: `research/drone_transmitters/06_dvbt_cofdm.md`, ETSI EN 300 744;
synthetic validation only, technology `cofdm_video_link` tagged `research`):
1. Welch spectrum candidates 0.951·B wide (B = 8 / 7 / 6 / 5 / 4 / 2.5 / 2 MHz).
2. Slice at 8B/7 MHz; CP fold for 2K / 8K × GI.
3. Integer offset from the continual pilots.
4. TPS: DBPSK, sync word, BCH(67,53), giving mode, GI, constellation, hierarchy, code rates, frame number and
   cell-ID byte.

`tables.hpp` is generated from python `tables.py`. Tests: `tests/unit/test_cofdm.cpp` against
`tests/vectors/cofdm/gen_cofdm.py` (synthetic DVB-T, no capture).
