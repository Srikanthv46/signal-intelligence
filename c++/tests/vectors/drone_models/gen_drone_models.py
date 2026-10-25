#!/usr/bin/env python3
"""Reference for core/drone_models (research, PLAN M0-M2): DJI product types 0-80, CTA-2063-A serial parsing, model
blocks and Remote ID model blocks.

    .venv/bin/python -I c++/tests/vectors/drone_models/gen_drone_models.py
"""
import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.core import drone_models as M   # noqa: E402

serials = ["1581F45T0000000000AB", "1596A1234567890", "1581F45T", "1581FO5T0000000000AB", "9999312A", "1748C123456789ABC",
           "1973512345", "", "15815"]
blocks = [["Walksnail Avatar HD", ["Mini 1S", "Nano V3", "Mini 1S"], ["tx_antennas 1"], "inferred", None],
          ["Walksnail Avatar HD", ["Moonlight"], [], "inferred", None], ["DJI", ["Mini 2"], ["x"], "decoded", None],
          ["DJI", ["Avata-like"], [], "inferred", None], ["X", [], [], "inferred", None], ["X", [], [], "decoded", "why"]]
rids = [{"uas_id": s, "ua_type": "helicopter or multirotor", **({"manufacturer": M.serial_manufacturer(s)}
                                                               if M.serial_manufacturer(s) else {})} for s in serials[:6]]
rids += [{"ua_type": 15}, {}]
serials += ["1581F6Z9A00000000001", "1581F446A00000000001", "1869AF11P1234", "1748FEV3HMQ12345", "1748CLTC00001"]
rids += [{"uas_id": s, "manufacturer": M.serial_manufacturer(s)} for s in ("1581F6Z9A00000000001", "1581F446A00000000001")]
ref = dict(dji=[[c, M.dji_product(c), M.dji_product_source(c)] for c in range(181)],
           serial_models=[[s, M.serial_model(s)] for s in serials], serials=[[s, M.serial_manufacturer(s)] for s in serials],
           blocks=[[b, M.model_block(*b)] for b in blocks], rids=[[r, M.remote_id_model(r)] for r in rids])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1))
print(len(ref["dji"]), len(ref["serials"]), len(ref["blocks"]), len(ref["rids"]))
