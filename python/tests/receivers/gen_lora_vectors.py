#!/usr/bin/env python3
"""Generate LoRa test vectors with gr-lora_sdr's transmitter (independent reference).

Run under radioconda's Python (gr-lora_sdr is not importable from .venv):
    ~/radioconda/bin/python tests/receivers/gen_lora_vectors.py

Writes testdata/lora/tx_<name>.cf32 (complex64, 250 kSps = 2x BW 125 kHz, packets
separated by silence) and tx_vectors.json (config + the payloads, in order).
"""
import json, os, time
import pmt
from gnuradio import gr, blocks, lora_sdr

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "testdata", "lora")
FS, BW = 250e3, 125e3
# (name, sf, cr, crc, payloads) -- LDRO follows the Semtech rule (mode 2): on for SF11/12 at 125 kHz
CONFIGS = [
    ("sf7_cr1", 7, 1, True, ["a", "hello", "0123456789abcdef0123456789abcdef", "x" * 60]),
    ("sf8_cr2_nocrc", 8, 2, False, ["no crc here", "ab"]),
    ("sf9_cr3", 9, 3, True, ["Hello LoRa 434", "18 byte payload!!!", "pkt#3"]),
    ("sf10_cr4", 10, 4, True, ["sf10 4/8", "z"]),
    ("sf11_cr1_ldro", 11, 1, True, ["ldro on"]),
]


class TX(gr.top_block):
    def __init__(self, sf, cr, crc, path):
        gr.top_block.__init__(self)
        self.tx = lora_sdr.lora_sdr_lora_tx(bw=int(BW), cr=cr, has_crc=crc, impl_head=False, samp_rate=int(FS),
                                            sf=sf, ldro_mode=2, frame_zero_padd=int(2 ** sf * FS / BW * 6))
        self.sink = blocks.file_sink(gr.sizeof_gr_complex, path)
        self.connect(self.tx, self.sink)


def main():
    os.makedirs(OUT, exist_ok=True)
    meta = []
    for name, sf, cr, crc, pays in CONFIGS:
        path = os.path.join(OUT, f"tx_{name}.cf32")
        tb = TX(sf, cr, crc, path)
        tb.start()
        for p in pays:
            tb.tx.lora_sdr_whitening_0.to_basic_block()._post(pmt.intern("msg"), pmt.intern(p))   # the hier "in" port is only an alias
            time.sleep(0.2)
        time.sleep(1.0)
        tb.stop(); tb.wait()
        meta.append(dict(name=name, sf=sf, bw=BW, fs=FS, cr=cr, crc=crc, payloads_hex=[p.encode().hex() for p in pays]))
    json.dump(meta, open(os.path.join(OUT, "tx_vectors.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
