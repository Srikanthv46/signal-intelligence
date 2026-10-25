# O4 transmit antennas (bench, PLAN D6 / M4b, 2026-10-22)

`probe.py`: `digital_vtx.streams.tx_streams` on the 13 O4 video bursts of `testdata/DJI_Mini_Pro_4_20251119-172214`
(Mini 4 Pro, 2T4R per DJI), slice at 30.72 MSps (Tu = 1024), 9.18 MHz, 295 symbols x 275 subcarriers:
kurtosis **1.70**, time-pair delta +0.00, frequency-pair delta -0.02 -> no verdict.

Reading: above any single-stream QAM (<= 1.40) and without the Alamouti pairing, i.e. spatial multiplexing,
DFT-spread data, or contamination by the O4 reference symbols (comb-2 references every 0.5 ms are in the runs).
Nothing here separates the O4 Air Unit (1 TX) from the Pro (2 TX): that needs the O4 reference / data symbol map
and a known single-antenna O4 capture. Not added to the library.
