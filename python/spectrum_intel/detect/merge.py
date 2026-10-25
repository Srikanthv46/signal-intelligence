"""Box clean-up: fuse the two resolutions, merge fragments, split adjacent
carriers at spectral valleys, and group bursts into per-emitter signals."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter1d

from ..model import Box, Signal
from .stft import Spectrogram

NARROW_HZ = 300e3       # below this the fine (3.75 kHz) spectrogram owns the box


def _overlap(a0, a1, b0, b1):
    return max(0.0, min(a1, b1) - max(a0, b0))


def fuse_resolutions(coarse: list[Box], fine: list[Box]) -> list[Box]:
    """Wide boxes from the coarse map; narrow ones from the fine map unless they
    lie inside a wide box. Narrow coarse boxes survive only where the fine map
    saw nothing (bursts shorter than its 267 us frames)."""
    wide = [b for b in coarse if b.bw >= NARROW_HZ]
    out = list(wide)
    fine_n = [b for b in fine if b.bw < NARROW_HZ]

    def inside_wide(b):
        return any(_overlap(b.f0, b.f1, w.f0, w.f1) > 0.5 * b.bw and
                   _overlap(b.t0, b.t1, w.t0, w.t1) > 0.5 * b.dur for w in wide)
    out += [b for b in fine_n if not inside_wide(b)]
    for b in coarse:
        if b.bw < NARROW_HZ and not inside_wide(b) and not any(
                _overlap(b.f0 - 10e3, b.f1 + 10e3, f.f0, f.f1) > 0 and
                _overlap(b.t0, b.t1, f.t0 - 300e-6, f.t1 + 300e-6) > 0 for f in fine_n):
            out.append(b)
    return out


def merge_fragments(boxes: list[Box], f_gap: float = 130e3, t_gap: float = 40e-6) -> list[Box]:
    """Union boxes that touch: same time and <= f_gap apart in frequency (a
    carrier split by its DC null / a deep fade), or same frequency span and
    <= t_gap apart in time (a burst broken by a fade)."""
    boxes = list(boxes)
    changed = True
    while changed:
        changed = False
        boxes.sort(key=lambda b: (b.f0, b.t0))
        out = []
        for b in boxes:
            for i, a in enumerate(out):
                tov = _overlap(a.t0, a.t1, b.t0, b.t1)
                fov = _overlap(a.f0, a.f1, b.f0, b.f1)
                side = (tov > 0.5 * min(a.dur, b.dur) and
                        max(a.f0, b.f0) - min(a.f1, b.f1) <= f_gap and a.res == b.res)
                seq = (fov > 0.8 * min(a.bw, b.bw) and max(a.t0, b.t0) - min(a.t1, b.t1) <= t_gap
                       and 0.6 < a.bw / b.bw < 1.6)
                if a.kind != b.kind:
                    continue
                if side or seq or (tov > 0 and fov > 0.5 * min(a.bw, b.bw)):
                    out[i] = Box(min(a.t0, b.t0), max(a.t1, b.t1), min(a.f0, b.f0), max(a.f1, b.f1),
                                 max(a.peak_snr_db, b.peak_snr_db),
                                 a.mean_snr_db if a.n_cells >= b.n_cells else b.mean_snr_db,
                                 a.res if a.bw >= b.bw else b.res, a.n_cells + b.n_cells, a.kind)
                    changed = True
                    break
            else:
                out.append(b)
        boxes = out
    return boxes


def split_valleys(b: Box, S: Spectrogram, floor_db: np.ndarray, depth_db: float = 6.0,
                  min_part_hz: float = 600e3) -> list[Box]:
    """Adjacent carriers (e.g. two operators' LTE blocks) touch in the energy
    map; cut a wide box where its mean spectrum dips >= depth_db below the
    plateaus on both sides."""
    if b.bw < 2 * min_part_hz:
        return [b]
    i0, i1 = S.f_index(b.f0 + S.df / 2), S.f_index(b.f1 - S.df / 2) + 1
    j0, j1 = S.t_index(b.t0), S.t_index(b.t1 - S.frame_len) + 1
    prof = 10 * np.log10(S.P[j0:j1, i0:i1].mean(0) + 1e-30) - floor_db[i0:i1]
    prof = uniform_filter1d(prof, 3)
    kmin = max(2, int(min_part_hz / S.df))
    cuts = []
    for k in range(kmin, len(prof) - kmin):
        if prof[k] != prof[max(0, k - 2):k + 3].min():
            continue
        left, right = prof[:k].max(), prof[k + 1:].max()
        lo = prof[max(0, k - 2 * kmin):k].max(), prof[k + 1:k + 1 + 2 * kmin].max()
        if min(lo) - prof[k] >= depth_db and min(left, right) - prof[k] >= depth_db:
            if not cuts or k - cuts[-1] >= kmin:
                cuts.append(k)
    if not cuts:
        return [b]
    edges = [0] + cuts + [len(prof)]
    out = []
    for a, c in zip(edges[:-1], edges[1:]):
        out.append(Box(b.t0, b.t1, float(S.f[i0 + a] - S.df / 2), float(S.f[i0 + c - 1] + S.df / 2),
                       round(float(prof[a:c].max()), 1), round(float(np.median(prof[a:c])), 1),
                       b.res, b.n_cells * (c - a) // len(prof), b.kind))
    return out


def group_signals(boxes: list[Box], iou_min: float = 0.6) -> list[Signal]:
    """Bursts whose frequency spans overlap (IoU >= iou_min) belong to one
    emitter/channel; the signal span is the median burst span."""
    boxes = sorted(boxes, key=lambda b: b.t0)
    groups: list[list[Box]] = []
    for b in boxes:
        best, best_iou = None, iou_min
        for g in groups:
            f0, f1 = np.median([x.f0 for x in g]), np.median([x.f1 for x in g])
            inter = _overlap(f0, f1, b.f0, b.f1)
            iou = inter / (max(f1, b.f1) - min(f0, b.f0))
            if iou >= best_iou:
                best, best_iou = g, iou
        if best is None:
            groups.append([b])
        else:
            best.append(b)
    sigs = []
    for g in groups:
        w = np.array([max(x.dur, 1e-6) for x in g])
        sigs.append(Signal(bursts=g, f0=float(np.median([x.f0 for x in g])),
                           f1=float(np.median([x.f1 for x in g])),
                           snr_db=round(float(10 * np.log10(np.average(
                               [10 ** (x.mean_snr_db / 10) for x in g], weights=w))), 1)))
    return sigs


def absorb_into_continuous(sigs: list[Signal], frac: float = 0.8) -> list[Signal]:
    """Load variation of an always-on carrier (e.g. LTE PDSCH) shows up as
    'bursts' above its 20th-percentile level; a burst-only signal whose band
    lies >= frac inside a continuous signal is part of that emitter."""
    cont = [s for s in sigs if any(b.kind == "continuous" for b in s.bursts)]
    out = []
    for s in sigs:
        if s in cont:
            out.append(s)
            continue
        host = next((c for c in cont if max(0.0, min(s.f1, c.f1) - max(s.f0, c.f0)) >= frac * s.bw), None)
        if host is None:
            out.append(s)
        else:
            host.bursts += s.bursts
            host.features["n_bursts_above_continuous_level"] = \
                host.features.get("n_bursts_above_continuous_level", 0) + len(s.bursts)
    return out


PT_PRI_TOL = 0.01        # pieces of one pulse train: same PRI within 1 %
PT_T_TOL_S = 5e-6        # and the same pulse starts within 5 us (1 us envelope + band-edge rise times)
PT_MATCH_FRAC = 0.5      # for at least half of the shorter train's pulses
PT_F_GAP_HZ = 2e6        # spans at most 2 MHz apart
PT_SLIVER_PAD_S = 20e-6  # a burst-only sliver whose every burst holds a pulse start (+-20 us) of the train


def _coincide(ta, tb) -> float:
    """Fraction of the shorter pulse list with a pulse of the other within PT_T_TOL_S."""
    a, b = (ta, tb) if len(ta) <= len(tb) else (tb, ta)
    if not a:
        return 0.0
    b = np.asarray(b)
    return sum(bool(np.any(np.abs(b - t) <= PT_T_TOL_S)) for t in a) / len(a)


def merge_pulse_trains(sigs: list[Signal]) -> tuple[list[Signal], list[Signal]]:
    """A wideband pulsed emitter (radar) leaves several narrower energy detections with one pulse train: the same
    PRI and the same pulse times. Merge them, and the burst-only slivers whose bursts all sit on its pulses, into
    one signal spanning their union. Returns (signals, merged): the merged signals need their features measured
    again."""
    hosts = [s for s in sigs if s.features.get("pri_s") and s.features.get("_pulse_t_s")]
    groups: list[list[Signal]] = []
    for s in hosts:
        g = next((g for g in groups if any(
            abs(s.features["pri_s"] - h.features["pri_s"]) <= PT_PRI_TOL * h.features["pri_s"]
            and max(s.f0, h.f0) - min(s.f1, h.f1) <= PT_F_GAP_HZ
            and _coincide(s.features["_pulse_t_s"], h.features["_pulse_t_s"]) >= PT_MATCH_FRAC for h in g)), None)
        if g is None:
            groups.append([s])
        else:
            g.append(s)
    used, merged = set(), []
    for g in groups:
        f0, f1 = min(h.f0 for h in g), max(h.f1 for h in g)
        t = np.asarray(sorted(v for h in g for v in h.features["_pulse_t_s"]))
        slivers = [s for s in sigs if all(s is not h for h in g) and id(s) not in used and not s.features.get("pri_s")
                   and not any(b.kind == "continuous" for b in s.bursts)
                   and s.f0 >= f0 - PT_F_GAP_HZ and s.f1 <= f1 + PT_F_GAP_HZ
                   and all(np.any((t >= b.t0 - PT_SLIVER_PAD_S) & (t <= b.t1 + PT_SLIVER_PAD_S)) for b in s.bursts)]
        if len(g) + len(slivers) < 2:
            continue
        parts = g + slivers
        top = max(g, key=lambda h: h.snr_db)
        m = Signal(bursts=[b for p in parts for b in p.bursts], f0=min(p.f0 for p in parts),
                   f1=max(p.f1 for p in parts), snr_db=top.snr_db, detection=dict(top.detection))
        m.detection["description"] = (m.detection.get("description", "") +
                                      f"; {len(parts)} energy detections with one pulse train (same PRI and pulse "
                                      "times) merged into one emitter")
        used.update(id(p) for p in parts)
        merged.append(m)
    out = [s for s in sigs if id(s) not in used] + merged
    return out, merged
