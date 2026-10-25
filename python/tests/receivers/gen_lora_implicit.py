"""Implicit-header LoRa vectors from gr-lora_sdr's transmitter (independent reference) for lora_phy.decode_implicit
(PLAN R2c): ExpressLRS sub-GHz shapes (BW 500 kHz, SF 6-9 (SX127x framing), CR 4/7 and 4/8, 8 and 13-byte payloads, LoRa CRC off).
ASCII payloads (gr-lora_sdr takes text).

Run under radioconda's Python:  ~/radioconda/bin/python tests/receivers/gen_lora_implicit.py
Writes testdata/lora/tx_implicit_<name>.cf32 (complex64, 1 MSps) and tx_implicit.json.
"""
import json, os, time
import pmt
from gnuradio import gr, blocks, lora_sdr

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "testdata", "lora")
FS, BW = 1e6, 500e3
CONFIGS = [
    ("sf6_cr3_len8", 6, 3, ["ELRS sf6", "200 Hz!!"]),
    ("sf7_cr3_len8", 7, 3, ["ELRSsync", "8 bytes!"]),
    ("sf8_cr3_len8", 8, 3, ["sf8 pkt1", "sf8 pkt2"]),
    ("sf9_cr3_len8", 9, 3, ["sf9 0123"]),
    ("sf7_cr4_len13", 7, 4, ["13 byte ota8!", "another 13 B."]),
]


class TX(gr.top_block):
    def __init__(self, sf, cr, n, path):
        gr.top_block.__init__(self)
        self.tx = lora_sdr.lora_sdr_lora_tx(bw=int(BW), cr=cr, has_crc=False, impl_head=True, samp_rate=int(FS),
                                            sf=sf, ldro_mode=2, frame_zero_padd=int(2 ** sf * FS / BW * 6))
        self.sink = blocks.file_sink(gr.sizeof_gr_complex, path)
        self.connect(self.tx, self.sink)


def main():
    os.makedirs(OUT, exist_ok=True)
    meta = []
    for name, sf, cr, pays in CONFIGS:
        path = os.path.join(OUT, f"tx_implicit_{name}.cf32")
        tb = TX(sf, cr, len(pays[0]), path)
        tb.start()
        for p in pays:
            tb.tx.lora_sdr_whitening_0.to_basic_block()._post(pmt.intern("msg"), pmt.intern(p))
            time.sleep(0.2)
        time.sleep(1.0)
        tb.stop(); tb.wait()
        meta.append(dict(name=name, sf=sf, bw=BW, fs=FS, cr=cr, crc=False, payloads_hex=[p.encode().hex() for p in pays]))
    json.dump(meta, open(os.path.join(OUT, "tx_implicit.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
