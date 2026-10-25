#!/usr/bin/env python3
"""C++ vs Python parity on testdata: same signals, labels, decode status and decoded fields.

    .venv/bin/python c++/tools/parity.py --cpp c++/build/si_analyze                 # everything
    .venv/bin/python c++/tools/parity.py --cpp ... --receivers lte,nr --inputs cellular
    .venv/bin/python c++/tools/parity.py --cpp ... --receivers none                  # energy pipeline only

Both implementations write the per-block report (python report.writer.build; C++
si_analyze --report) for the same input and receivers. Signals are matched by label and
decode status, then by centre frequency and time; every decoded field, the label,
confidence, modulation, SNR, span, flags and the report-level extras (emitters, devices,
findings) are compared. Numbers may differ by one unit in the last decimal place the
Python value shows (rounded report values), or by 1e-6 relative.

Python reports are cached in out/parity_cache/<receivers>/ (delete it after changing the
Python code). Exit status 1 if any input differs.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TESTDATA = REPO / "testdata"
# receiver frequency error given with an input group (the Pluto's TCXO, measured from the LTE PSS at 940.1 MHz)
PPM = {"cellular": 7.44, "live": 7.4}
CACHE = REPO / "out" / "parity_cache"
IGNORE_TOP = {"source_path", "processing", "file"}


# ----------------------------------------------------------------------------- inputs
def inputs() -> list[dict]:
    """Every testdata input: SigMF captures, plus raw .cf32 vectors with a .json sidecar."""
    out = []
    for m in sorted(TESTDATA.rglob("*.sigmf-meta")):
        grp = m.parent.parent.name if m.parent.parent != TESTDATA else m.parent.name
        # a capture's own measured receiver error (bench:receiver_ppm) overrides its group's
        own = json.loads(m.read_text())["global"].get("bench:receiver_ppm")
        out.append({"name": m.stem.replace(".sigmf-meta", ""), "kind": "sigmf", "path": str(m), "group": grp,
                    "ppm": own if own is not None else PPM.get(grp)})
    for j in sorted(TESTDATA.rglob("*.json")):
        d = j.with_suffix(".cf32")
        if not d.exists():
            continue
        meta = json.loads(j.read_text())
        fs = meta.get("fs") or meta.get("fs_hz") or meta.get("sample_rate_hz")
        fc = meta.get("fc") or meta.get("fc_hz") or meta.get("center_freq_hz")
        if not fs or fc is None:
            continue
        out.append({"name": d.stem, "kind": "raw", "path": str(d), "fs": float(fs), "fc": float(fc),
                    "group": j.parent.name, "ppm": PPM.get(j.parent.name)})
    return out


# ----------------------------------------------------------------------------- runners
def python_report(inp: dict, receivers: list[str] | None) -> dict:
    import numpy as np
    from spectrum_intel.analyze import analyze
    from spectrum_intel.decode import registry
    from spectrum_intel.io.sigmf import Capture
    from spectrum_intel.report.writer import build
    plugins = [] if receivers == [] else registry.plugins(receivers)
    ppm = {"ppm": inp["ppm"], "source": "caller (IQBlock.receiver_ppm)"} if inp.get("ppm") is not None else None
    if inp["kind"] == "sigmf":
        a = analyze(inp["path"], None, plugins, ppm=ppm)
    else:
        x = np.fromfile(inp["path"], np.complex64)
        a = analyze(Capture.from_iq(x, inp["fs"], inp["fc"], name=inp["name"]), None, plugins, ppm=ppm)
    for s in a.signals:
        s.features.pop("_video_frame", None)
        s.features.pop("_audio", None)
    return json.loads(json.dumps(build(a), default=str))


def cpp_report(cpp: str, inp: dict, receivers: list[str] | None, threads: int) -> tuple[dict, float]:
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "r.json"
        cmd = [cpp, "--report", "--out", str(out), "--indent", "-1"]
        if threads:
            cmd += ["--threads", str(threads)]
        if inp.get("ppm") is not None:
            cmd += ["--ppm", repr(inp["ppm"])]
        if receivers is not None:
            cmd += ["--enable", ",".join(receivers) if receivers else "none"]
        if inp["kind"] == "sigmf":
            cmd.append(inp["path"])
        else:
            cmd += ["--raw", inp["path"], "--format", "cf32", "--fs", repr(inp["fs"]), "--fc", repr(inp["fc"])]
        t = time.time()
        p = subprocess.run(cmd, capture_output=True, text=True)
        dt = time.time() - t
        if p.returncode != 0 or not out.exists():
            raise RuntimeError(f"si_analyze failed ({p.returncode}): {p.stderr[-2000:]}")
        return json.loads(out.read_text()), dt


# ----------------------------------------------------------------------------- comparison
def _decimals(x) -> int:
    s = repr(x)
    if "e" in s or "E" in s:
        return 15
    return len(s.split(".")[1]) if "." in s else 0


def _num_ok(a, b) -> bool:
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b
    if isinstance(a, float) and math.isnan(a):
        return isinstance(b, float) and math.isnan(b)
    tol = max(1e-6 * abs(a), 10.0 ** -_decimals(a) if isinstance(a, float) else 0.0)
    return abs(a - b) <= tol * (1 + 1e-9)


def diff(a, b, path: str, out: list[str]):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a.keys() | b.keys():
            if k not in a:
                out.append(f"{path}.{k}: only in C++ ({str(b[k])[:60]})")
            elif k not in b:
                out.append(f"{path}.{k}: missing in C++ (python {str(a[k])[:60]})")
            else:
                diff(a[k], b[k], f"{path}.{k}", out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: {len(a)} items in python, {len(b)} in C++")
        for i, (x, y) in enumerate(zip(a, b)):
            diff(x, y, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not _num_ok(a, b):
            out.append(f"{path}: python {a!r} != C++ {b!r}")
    elif a != b:
        out.append(f"{path}: python {str(a)[:80]!r} != C++ {str(b)[:80]!r}")


def _key(s: dict) -> str:
    return (s.get("features") or {}).get("decoder") or (s["decoded"].get("method") or "energy").split(" ")[0]


def match(py: list[dict], cc: list[dict]) -> tuple[list, list, list]:
    """Pairs (python, C++) by label + status, then nearest centre / time overlap."""
    pairs, used = [], set()
    for s in py:
        best, bd = None, None
        for j, c in enumerate(cc):
            if j in used or c["label"] != s["label"] or c["decoded"]["status"] != s["decoded"]["status"]:
                continue
            tov = min(s["t_end_s"], c["t_end_s"]) - max(s["t_start_s"], c["t_start_s"])
            df = abs(s["f_center_hz"] - c["f_center_hz"])
            if df > max(50e3, 0.05 * max(s["bandwidth_hz"], 1)) or tov < -1e-4:
                continue
            if bd is None or df < bd:
                best, bd = j, df
        if best is None:
            pairs.append((s, None))
        else:
            used.add(best)
            pairs.append((s, cc[best]))
    extra = [c for j, c in enumerate(cc) if j not in used]
    return [p for p in pairs if p[1] is not None], [p[0] for p in pairs if p[1] is None], extra


def compare(py: dict, cc: dict) -> dict:
    res = {"n_python": len(py["signals"]), "n_cpp": len(cc["signals"]), "diffs": [], "missing": [], "extra": [],
           "by_receiver": {}}
    pairs, missing, extra = match(py["signals"], cc["signals"])
    for s, c in pairs:
        d = []
        a = {k: v for k, v in s.items() if k != "id"}
        b = {k: v for k, v in c.items() if k != "id"}
        diff(a, b, f"signal[{s['label'][:40]} @ {s['f_center_hz'] / 1e6:.3f} MHz]", d)
        res["diffs"] += d
        r = res["by_receiver"].setdefault(_key(s), {"matched": 0, "identical": 0, "missing": 0, "extra": 0})
        r["matched"] += 1
        r["identical"] += not d
    for s in missing:
        res["missing"].append(f"{s['label']} @ {s['f_center_hz'] / 1e6:.4f} MHz [{s['decoded']['status']}]")
        res["by_receiver"].setdefault(_key(s), {"matched": 0, "identical": 0, "missing": 0, "extra": 0})["missing"] += 1
    for c in extra:
        res["extra"].append(f"{c['label']} @ {c['f_center_hz'] / 1e6:.4f} MHz [{c['decoded']['status']}]")
        res["by_receiver"].setdefault(_key(c), {"matched": 0, "identical": 0, "missing": 0, "extra": 0})["extra"] += 1
    top = []
    diff({k: v for k, v in py.items() if k not in IGNORE_TOP | {"signals"}},
         {k: v for k, v in cc.items() if k not in IGNORE_TOP | {"signals"}}, "report", top)
    res["diffs"] += top
    res["ok"] = not (res["diffs"] or res["missing"] or res["extra"])
    return res


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cpp", default=str(REPO / "c++" / "build" / "si_analyze"), help="si_analyze binary")
    ap.add_argument("--receivers", default=None, help="comma list for both sides; 'none' = energy pipeline only")
    ap.add_argument("--inputs", default=None, help="comma list of input names or groups (default: all)")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--show", type=int, default=12, help="differences shown per input")
    ap.add_argument("--refresh", action="store_true", help="recompute the cached Python reports")
    ap.add_argument("--save", default=None, help="write the C++ reports to this folder")
    ap.add_argument("--list", action="store_true", help="list the inputs and exit")
    ap.add_argument("--python-only", action="store_true", help="only fill the Python report cache")
    args = ap.parse_args()

    rx = None if args.receivers is None else ([] if args.receivers == "none" else args.receivers.split(","))
    sel = inputs()
    if args.list:
        for i in sel:
            print(f"{i['group']:>40} / {i['name']} ({i['kind']})")
        return 0
    if args.inputs:
        want = set(args.inputs.split(","))
        sel = [i for i in sel if i["name"] in want or i["group"] in want]
    tag = "all" if rx is None else ("none" if rx == [] else "-".join(sorted(rx)))
    cache = CACHE / tag
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("SPECTRUM_INTEL_THREADS", str(args.threads))

    bad = 0
    totals: dict[str, dict] = {}
    for inp in sel:
        cp = cache / f"{inp['name']}.json"
        t = time.time()
        if cp.exists() and not args.refresh:
            py = json.loads(cp.read_text())
            tpy = None
        else:
            py = python_report(inp, rx)
            tpy = time.time() - t
            cp.write_text(json.dumps(py))
        if args.python_only:
            print(f"cached {inp['name']}" + (f" ({tpy:.1f} s)" if tpy else ""))
            continue
        try:
            cc, tcc = cpp_report(args.cpp, inp, rx, args.threads)
        except Exception as e:                                        # noqa: BLE001
            print(f"FAIL {inp['name']}: {e}")
            bad += 1
            continue
        if args.save:
            Path(args.save).mkdir(parents=True, exist_ok=True)
            (Path(args.save) / f"{inp['name']}_report.json").write_text(json.dumps(cc, indent=1))
        r = compare(py, cc)
        bad += not r["ok"]
        speed = f"C++ {tcc:.2f} s" + (f", python {tpy:.1f} s" if tpy else "")
        print(f"{'OK  ' if r['ok'] else 'DIFF'} {inp['name']}: python {r['n_python']} / C++ {r['n_cpp']} signals; "
              f"{len(r['diffs'])} field diffs, {len(r['missing'])} missing, {len(r['extra'])} extra ({speed})")
        for line in (r["diffs"][:args.show] + [f"missing in C++: {m}" for m in r["missing"][:args.show]]
                     + [f"extra in C++: {m}" for m in r["extra"][:args.show]]):
            print("     ", line)
        for k, v in r["by_receiver"].items():
            t_ = totals.setdefault(k, {"matched": 0, "identical": 0, "missing": 0, "extra": 0})
            for kk in t_:
                t_[kk] += v[kk]
    print("\nper receiver (matched / identical / missing in C++ / extra in C++):")
    for k, v in sorted(totals.items()):
        print(f"  {k:>18}: {v['matched']:4d} / {v['identical']:4d} / {v['missing']:3d} / {v['extra']:3d}")
    print(f"\n{len(sel) - bad} of {len(sel)} inputs identical within tolerance")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
