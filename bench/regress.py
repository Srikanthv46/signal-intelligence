#!/usr/bin/env python3
"""Bench: Python regression harness over testdata (all receivers, the same runner as c++/tools/parity.py).

    .venv/bin/python bench/regress.py run  <tag> [--inputs a,b] [--workers 2]   -> out/regress/<tag>/<input>.json
    .venv/bin/python bench/regress.py diff <tagA> <tagB>                         -> per-input signal changes
"""
import argparse, json, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "c++" / "tools"))
import parity  # noqa: E402

OUT = REPO / "out" / "regress"


def _inputs(sel):
    ins = parity.inputs()                     # receiver ppm per group / per capture (parity.PPM, bench:receiver_ppm)
    if sel:
        want = set(sel.split(","))
        ins = [i for i in ins if i["name"] in want or i["group"] in want]
    return ins


def _one(args):
    inp, d = args
    t = time.time()
    r = parity.python_report(inp, None)
    r["_regress_s"] = round(time.time() - t, 2)
    (d / f"{inp['name']}.json").write_text(json.dumps(r))
    return inp["name"], r["_regress_s"], len(r.get("signals", []))


def run(tag, sel, workers):
    d = OUT / tag; d.mkdir(parents=True, exist_ok=True)
    ins = _inputs(sel)
    with ProcessPoolExecutor(workers) as ex:
        for name, dt, n in ex.map(_one, [(i, d) for i in ins]):
            print(f"{name:45s} {dt:7.1f} s  {n:3d} signals", flush=True)


def _sig(s):
    f = s.get("f_center_hz") or s.get("frequency", {}).get("center_hz")
    d = s.get("decoded") or s.get("decode") or {}
    return (round((f or 0) / 1e3), str(s.get("label")), str(d.get("status")),
            ",".join(sorted(map(str, s.get("flags") or []))))


def diff(a, b):
    A, B = OUT / a, OUT / b
    for fa in sorted(A.glob("*.json")):
        fb = B / fa.name
        if not fb.exists():
            continue
        ra, rb = json.loads(fa.read_text()), json.loads(fb.read_text())
        sa, sb = {_sig(s) for s in ra["signals"]}, {_sig(s) for s in rb["signals"]}
        gone, new = sorted(sa - sb, key=lambda x: x[0]), sorted(sb - sa, key=lambda x: x[0])
        print(f"{fa.stem:45s} {ra.get('_regress_s', 0):6.1f}s -> {rb.get('_regress_s', 0):6.1f}s  "
              f"signals {len(ra['signals'])} -> {len(rb['signals'])}  changed {len(gone)}/{len(new)}")
        for g in gone:
            print(f"     - {g[0]/1e3:11.3f} MHz  {g[1]}  [{g[2]}] {g[3] or ''}")
        for n in new:
            print(f"     + {n[0]/1e3:11.3f} MHz  {n[1]}  [{n[2]}] {n[3] or ''}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("tags", nargs="+")
    ap.add_argument("--inputs"); ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    run(a.tags[0], a.inputs, a.workers) if a.cmd == "run" else diff(*a.tags)
