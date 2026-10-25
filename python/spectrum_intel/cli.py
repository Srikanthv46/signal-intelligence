"""spectrum-intel: detect, identify and decode every signal in SigMF snapshots.

    spectrum-intel run PATH [--out DIR] [--workers N] [--enable R,..] [--disable R,..]
    spectrum-intel receivers                 list installed receivers

PATH is a snapshot (.sigmf-meta / its folder) or a sweep folder. With a sweep
folder every snapshot is processed and a sweep summary is written. A single
snapshot still uses its sweep siblings (if any) for the noise template and
the receiver ppm calibration.

Outputs (OUT mirrors snapshot names):
    <name>_report.json, <name>_annotated.png        per snapshot
    summary.json, signals.csv, occupancy_map.png    per sweep
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .io.sigmf import find_captures, load

DEFAULT_OUT = Path(__file__).resolve().parents[2] / "out"                # <repo>/out


def _sweep_root(p: Path) -> Path | None:
    """Folder whose sub-folders hold sibling snapshots of the same sweep (signal-vault
    names them SWEEP_*; a folder of unrelated captures is not a sweep)."""
    for cand in (p.parent.parent, p.parent):
        if _is_sweep(cand, list(cand.glob("*/*.sigmf-meta")) + list(cand.glob("*.sigmf-meta"))):
            return cand
    return None


def _is_sweep(folder: Path, metas: list[Path]) -> bool:
    return len(metas) > 4 and ("SWEEP" in folder.name or sum("SWEEP" in m.name for m in metas) > len(metas) / 2)


_SWEEP = None
_PLUGINS = None


def _init(sweep_cache, enable=None, disable=None):
    global _SWEEP, _PLUGINS
    from .decode import registry
    from .sweep import SweepContext
    _SWEEP = SweepContext.load(Path(sweep_cache)) if sweep_cache else None
    _PLUGINS = registry.plugins(enable, disable)


def process_one(meta: str, out_dir: str, png: bool = True) -> dict:
    from .analyze import analyze
    from .report import writer
    from .viz import annotate
    t = time.time()
    try:
        a = analyze(meta, _SWEEP, _PLUGINS)
        od = Path(out_dir)
        od.mkdir(parents=True, exist_ok=True)
        _write_pictures(a, od)
        _write_audio(a, od)
        rep = writer.build(a)
        rp = od / f"{a.cap.name}_report.json"
        rp.write_text(json.dumps(rep, indent=1))
        if png:
            annotate.draw(rep, a.cap.iq, od / f"{a.cap.name}_annotated.png")
        return {"meta": meta, "report": str(rp), "n": len(rep["signals"]), "s": round(time.time() - t, 1)}
    except Exception as e:                                  # noqa: BLE001 — one bad file must not stop a sweep
        return {"meta": meta, "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-1500:]}


def _write_pictures(a, od: Path):
    """Demodulated analog video -> <name>_vtx<MHz>_video.png, path recorded in decoded.fields."""
    from PIL import Image
    for s in a.signals:
        img = s.features.pop("_video_frame", None)
        if img is None:
            continue
        f_mhz = s.decoded["fields"].get("center_mhz") or s.fc / 1e6
        fn = f"{a.cap.name}_vtx{f_mhz:.1f}MHz_video.png"
        Image.fromarray(img).resize((img.shape[1], img.shape[0])).save(od / fn)
        s.decoded["fields"]["video_frame_png"] = fn


def _write_audio(a, od: Path):
    """Demodulated FM voice -> <name>_voice<MHz>_audio.wav (16-bit mono), path recorded in decoded.fields."""
    import wave
    import numpy as np
    for s in a.signals:
        au = s.features.pop("_audio", None)
        if au is None:
            continue
        x, fs = au
        f_mhz = s.decoded["fields"].get("center_mhz") or s.fc / 1e6
        fn = f"{a.cap.name}_voice{f_mhz:.4f}MHz_audio.wav"
        with wave.open(str(od / fn), "wb") as w:
            w.setnchannels(1)
            w.setsampwidth(2)
            w.setframerate(int(round(fs)))
            w.writeframes((np.clip(x, -1, 1) * 32767).astype("<i2").tobytes())
        s.decoded["fields"]["audio_wav"] = fn


def run(path: Path, out: Path, workers: int, png: bool = True, summary: bool = True, resume: bool = False,
        freq_range: list[tuple[float, float]] | None = None, enable: list[str] | None = None,
        disable: list[str] | None = None):
    from .decode import registry
    from .sweep import build_context
    path = path.resolve()
    sel = registry.plugins(enable, disable)                 # validates the names before any work
    print(f"plugins: {', '.join(p.name for p in sel)}"
          + (f"  (unavailable: {registry.unavailable()})" if registry.unavailable() else ""), flush=True)
    if path.is_dir() and not list(path.glob("*.sigmf-meta")):
        metas = find_captures(path)
        root = path if _is_sweep(path, metas) else None
        sweep_out = out / path.name
    else:
        metas, root = [load(path).path], _sweep_root(load(path).path)
        sweep_out = out / (root.name if root else "single")
    if freq_range:
        metas = [m for m in metas if any(lo <= load(m).fc <= hi for lo, hi in freq_range)]
        print(f"freq range {', '.join(f'{lo / 1e6:g}-{hi / 1e6:g}' for lo, hi in freq_range)} MHz "
              f"(labelled centre): {len(metas)} snapshot(s)", flush=True)
    cache = None
    if root:
        print(f"sweep context from {root} ...", flush=True)
        ctx = build_context(root, sweep_out / "_sweep", workers)
        cache = str(sweep_out / "_sweep" / "sweep_context.npz")
        if ctx.ppm:
            print(f"  receiver ppm {ctx.ppm.get('ppm')}", flush=True)
    if resume:
        done = {p.name[:-len("_report.json")] for p in sweep_out.rglob("*_report.json")}
        skipped = [m for m in metas if m.name[:-len(".sigmf-meta")] in done]
    else:
        skipped = []
    if resume:
        metas = [m for m in metas if m.name[:-len(".sigmf-meta")] not in done]
        print(f"resume: {len(skipped)} snapshot(s) already have reports", flush=True)
    print(f"processing {len(metas)} snapshot(s) -> {sweep_out}", flush=True)
    results = []
    t0 = time.time()
    if len(metas) == 1 or workers <= 1:
        _init(cache, enable, disable)
        for m in metas:
            results.append(process_one(str(m), str(sweep_out / m.parent.name), png))
            _log(results[-1], len(results), len(metas))
    else:
        with ProcessPoolExecutor(workers, initializer=_init, initargs=(cache, enable, disable)) as ex:
            futs = [ex.submit(process_one, str(m), str(sweep_out / m.parent.name), png) for m in metas]
            for f in as_completed(futs):
                results.append(f.result())
                _log(results[-1], len(results), len(metas))
    print(f"done in {time.time() - t0:.0f} s, {sum('error' in r for r in results)} failed", flush=True)
    log = sweep_out / "run_log.json"
    if resume:
        # every report on disk, whichever run produced it
        have = {r["meta"] for r in results}
        for m in skipped:
            if str(m) not in have:
                rp = sweep_out / m.parent.name / (m.name[:-len(".sigmf-meta")] + "_report.json")
                results.append({"meta": str(m), "report": str(rp)})
    if log.exists():                                        # merge: a partial rerun must not drop the others
        try:
            old = {r["meta"]: r for r in json.loads(log.read_text())}
        except (ValueError, KeyError, TypeError):
            old = {}
        old.update({r["meta"]: r for r in results})
        results_all = list(old.values())
    else:
        results_all = results
    log.write_text(json.dumps(results_all, indent=1))
    if summary and len(results) > 1:
        from . import summary as S
        S.build(sweep_out, png=png)
    return results


def _log(r, i, n):
    if "error" in r:
        print(f"  [{i}/{n}] FAILED {Path(r['meta']).name}: {r['error']}", flush=True)
    else:
        print(f"  [{i}/{n}] {r['n']:3d} signals {r['s']:5.1f}s {Path(r['meta']).name}", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="spectrum-intel", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="analyse a snapshot or a sweep folder")
    r.add_argument("path", type=Path)
    r.add_argument("--out", type=Path, default=DEFAULT_OUT)
    r.add_argument("--workers", type=int, default=2)        # 8 GB RAM: ~1 GB per worker
    r.add_argument("--no-png", action="store_true")
    r.add_argument("--resume", action="store_true", help="skip snapshots that already have a report")
    r.add_argument("--freq-range", metavar="LO:HI", help="only snapshots whose labelled centre is in LO..HI MHz (comma-separate several)")
    r.add_argument("--enable", metavar="R1,R2", help="run only these receivers (see `spectrum-intel receivers`)")
    r.add_argument("--disable", metavar="R1,R2", help="skip these receivers")
    sub.add_parser("receivers", help="list the receivers (plugins) that are installed")
    s = sub.add_parser("summary", help="(re)build the sweep summary from existing reports")
    s.add_argument("sweep_out", type=Path)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        fr = [tuple(float(v) * 1e6 for v in r.split(":")) for r in a.freq_range.split(",")] if a.freq_range else None
        csv = lambda v: [x.strip() for x in v.split(",") if x.strip()] if v else None   # noqa: E731
        run(a.path, a.out, a.workers, png=not a.no_png, resume=a.resume, freq_range=fr,
            enable=csv(a.enable), disable=csv(a.disable))
    elif a.cmd == "receivers":
        import importlib
        from .decode import registry
        for n in registry.receiver_names():
            doc = (importlib.import_module(f"spectrum_intel.receivers.{n}").__doc__ or "").strip().splitlines()
            names = [p.name for p in registry.plugins(enable=[n])]
            print(f"{n:12s} plugins {','.join(names) or '-':22s} {doc[0] if doc else ''}")
        for n, err in registry.unavailable().items():
            print(f"{n:12s} UNAVAILABLE: {err}")
    else:
        from . import summary as S
        S.build(a.sweep_out)


if __name__ == "__main__":
    sys.exit(main())
