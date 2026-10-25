"""Pulse deinterleaving (research, signal gaps S10; research/signal_gaps/05_aviation_radar_jammers.md sections
4.2-4.3): split one signal's pulse list (TOA, PW) into emitters by PRI. Deterministic; pure Python on lists;
mirrored in C++ core/deinterleave.

    deinterleave(toa_s, pw_s) -> {"emitters": [...], "unassigned": n}

1. PW groups: pulses sorted by width; a group holds widths within max(1.5 us, 30 %) of its narrowest.
2. Per group, up to 4 emitters, each:
   - CDIF: the differences t[i+c] - t[i] of levels 1..c (c = 1..8, accumulated) in 20 us - 20 ms, sorted and
     clustered (consecutive within max(1.5 us, 2 %)); a cluster with >= max(3, 0.25 (E - c)) members (E pulses left)
     is a PRI candidate (its median, with the level it appeared at).
   - Sequence search per candidate, strict (tolerance max(1.5 us, 2 % / sqrt(c)): a frame period sums c intervals)
     and in jitter mode (10 %, one chain; kept only when it explains > 1.25 x the pulses of the best strict
     emitter, since a stagger also fits a 10 % window): from
     every pulse not yet in an explored chain, step by the PRI with up to 2 missed pulses; the longest chain is kept
     when >= 60 % of the pulses its span implies are in it and random pulses at this density would form it with
     p < 1e-4 per block (the product over its steps of k x rate x 2 tol, k = 1 + missed pulses at that step, times
     E x 8 starts / candidates). Chains are taken until none is left; those overlapping the first in time are one
     emitter (>= 2 of them: a stagger, frame period = the PRI, sub-intervals from their phases).
   - The candidate whose emitter explains the most pulses (minus missed slots, minus 6 per extra chain) wins; its
     pulses are removed. A later strict emitter at the same period joins it when it overlaps in time (more
     stagger phases) or sits on its PRI grid (the same train after a beam-scan gap).
   - Phases (chain starts modulo the period) within the tolerance are one; a "stagger" whose sub-intervals are all
     equal is a constant PRI (the frame / K).
3. An emitter with >= 6 pulses is reported: pri_us (the frame period for a stagger), n_pulses, pw_us, stagger_us
   (sub-intervals in order) or None, jitter (rms of the in-chain intervals around the PRI / PRI), t_first_s.
   The rest stay unassigned.
"""
from __future__ import annotations

import math

MIN_D, MAX_D = 20e-6, 20e-3
TOL_ABS, TOL_REL = 1.5e-6, 0.02
PW_ABS, PW_REL = 1.5e-6, 0.30
MAX_LEVEL = 8
CDIF_X = 0.25
MAX_MISS = 2
MIN_CHAIN = 4
MIN_FILL = 0.6
MAX_EMITTERS = 4
MAX_CAND = 12
P_CHANCE = 1e-4
MIN_REPORT = 6
CHAIN_COST = 6.0
TOL_JIT = 0.10
JIT_GAIN = 1.25


def _median(v: list[float]) -> float:
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def cdif(r: list[float]) -> list[tuple[float, int]]:
    """(PRI candidate, level) of the sorted TOAs r, by level, then ascending; at most MAX_CAND."""
    e = len(r)
    out: list[tuple[float, int]] = []
    acc: list[float] = []
    for c in range(1, min(MAX_LEVEL, e - 1) + 1):
        acc += [r[i + c] - r[i] for i in range(e - c) if MIN_D <= r[i + c] - r[i] <= MAX_D]
        d = sorted(acc)
        need = max(3.0, CDIF_X * (e - c))
        i = 0
        while i < len(d):
            j = i + 1
            while j < len(d) and d[j] - d[j - 1] <= max(TOL_ABS, TOL_REL * d[j - 1]):
                j += 1
            if j - i >= need:
                m = _median(d[i:j])
                if all(abs(m - q) > max(TOL_ABS, TOL_REL * q) for q, _ in out):
                    out.append((m, c))
            i = j
        if len(out) >= MAX_CAND:
            break
    return out[:MAX_CAND]


def chain(r: list[float], pri: float, tol: float, rate: float, n_all: int) -> tuple[list[int], int]:
    """(indices into r, missed slots) of the longest chain at this PRI that is not chance, or ([], 0)."""
    best: list[int] = []
    best_k: list[int] = []
    seen = [False] * len(r)
    for s in range(len(r)):
        if seen[s]:
            continue
        cur, ks = [s], []
        j = s
        while True:
            nxt, kk = None, 0
            for k in range(1, MAX_MISS + 2):
                want = r[cur[-1]] + k * pri
                q = j + 1
                while q < len(r) and r[q] < want - tol:
                    q += 1
                if q < len(r) and abs(r[q] - want) <= tol:
                    nxt, kk = q, k
                    break
            if nxt is None:
                break
            cur.append(nxt)
            ks.append(kk)
            seen[nxt] = True
            j = nxt
        if len(cur) > len(best):
            best, best_k = cur, ks
    if len(best) < MIN_CHAIN:
        return [], 0
    span = r[best[-1]] - r[best[0]]
    if len(best) < MIN_FILL * (span / pri + 1):
        return [], 0
    logp = math.log(n_all * MAX_LEVEL) + sum(math.log(min(1.0, max(k * rate * 2 * tol, 1e-300))) for k in best_k)
    if logp >= math.log(P_CHANCE):
        return [], 0
    return best, sum(k - 1 for k in best_k)


def _try(r: list[float], w: list[float], pri: float, tol: float, single: bool = False):
    """The emitter (chains overlapping the first; only the first when single) this candidate explains at this
    tolerance: (score, chains, raw = pulses minus missed slots)."""
    rate = (len(r) - 1) / (r[-1] - r[0]) if r[-1] > r[0] else 0.0
    chains = []
    rem_t, rem_w = list(r), list(w)
    while len(rem_t) >= MIN_CHAIN and not (single and chains):
        c, miss = chain(rem_t, pri, tol, rate, len(r))
        if not c:
            break
        chains.append(([rem_t[i] for i in c], [rem_w[i] for i in c], miss))
        keep = set(c)
        rem_t = [t for i, t in enumerate(rem_t) if i not in keep]
        rem_w = [x for i, x in enumerate(rem_w) if i not in keep]
    if not chains:
        return None
    lo, hi = chains[0][0][0], chains[0][0][-1]
    used = [ch for ch in chains if ch[0][0] < hi and ch[0][-1] > lo]
    raw = sum(len(ch[0]) - ch[2] for ch in used)
    return raw - CHAIN_COST * (len(used) - 1), used, raw


def _same_phase(a: float, b: float, pri: float) -> bool:
    """b continues the train through a on the PRI grid (a beam-scan gap broke the chain)."""
    n = round((b - a) / pri)
    return abs(b - a - n * pri) <= max(TOL_ABS, TOL_REL * pri)


def _describe(pri: float, used: list) -> dict:
    ts = sorted(t for ch in used for t in ch[0])
    ws = [x for ch in used for x in ch[1]]
    lo = min(ch[0][0] for ch in used)
    tol = max(TOL_ABS, TOL_REL * pri / len(used))
    ph: list[float] = []                         # distinct phases: a continuation chain repeats one
    for x in sorted((ch[0][0] - lo) % pri for ch in used):
        if not ph or x - ph[-1] > tol:
            ph.append(x)
    if len(ph) >= 2 and ph[0] + pri - ph[-1] <= tol:
        ph.pop()                                 # the circular wrap
    K = len(ph)
    sub = [b - a for a, b in zip(ph, ph[1:] + [ph[0] + pri])]
    stagger = None
    if K >= 2 and max(sub) - min(sub) > 2 * tol:
        stagger = [round(x * 1e6, 2) for x in sub]
    elif K >= 2:
        pri = pri / K                            # equal steps: a constant PRI seen as a frame of K chains
    q = pri * K if stagger is None and K >= 2 else pri
    steps = []
    for ch in used:
        steps += [(b - a) / round((b - a) / q) for a, b in zip(ch[0][:-1], ch[0][1:])]
    jit = math.sqrt(sum((x - q) ** 2 for x in steps) / len(steps)) / q if steps else 0.0
    return {"pri_us": round(pri * 1e6, 2), "n_pulses": len(ts), "pw_us": round(_median(ws) * 1e6, 2),
            "stagger_us": stagger, "jitter": round(jit, 4), "t_first_s": round(ts[0], 7)}


def _emitters(r: list[float], w: list[float]) -> tuple[list[dict], list[float]]:
    found: list[tuple[float, list, str]] = []
    while len(r) >= MIN_CHAIN and len(found) < MAX_EMITTERS:
        best = jit = None
        cand = cdif(r)
        for pri, level in cand:                  # strict: a frame period sums `level` intervals
            t = _try(r, w, pri, max(TOL_ABS, TOL_REL * pri / math.sqrt(level)))
            if t is not None and (best is None or t[0] > best[0]):
                best = (t[0], t[1], pri, t[2])
        for pri, _ in cand:                      # jitter mode: one chain at 10 %
            t = _try(r, w, pri, max(TOL_ABS, TOL_JIT * pri), single=True)
            if t is not None and (jit is None or t[0] > jit[0]):
                jit = (t[0], t[1], pri, t[2])
        # a stagger also fits a 10 % window: jitter only when it explains clearly more than any strict emitter
        mode = "strict"
        if jit is not None and (best is None or jit[3] > JIT_GAIN * best[3]):
            best, mode = jit, "jitter"
        if best is None:
            break
        _, used, pri, _ = best
        prev = next((f for f in found if mode == "strict" and f[2] == "strict"
                     and abs(f[0] - pri) <= max(TOL_ABS, TOL_REL * pri)
                     and (min(ch[0][0] for ch in f[1]) < max(ch[0][-1] for ch in used)
                          and max(ch[0][-1] for ch in f[1]) > min(ch[0][0] for ch in used)
                          or _same_phase(f[1][0][0][0], used[0][0][0], pri))), None)
        if prev is not None:                     # more phases of a stagger, or the same train after a gap
            prev[1].extend(used)
        else:
            found.append((pri, used, mode))
        gone = {t for ch in used for t in ch[0]}
        keep = [i for i, t in enumerate(r) if t not in gone]
        r, w = [r[i] for i in keep], [w[i] for i in keep]
    return [_describe(pri, used) for pri, used, _ in found], r


def deinterleave(toa_s: list[float], pw_s: list[float]) -> dict:
    order = sorted(range(len(toa_s)), key=lambda i: (pw_s[i], toa_s[i]))
    groups: list[list[int]] = []
    for i in order:
        if groups and pw_s[i] <= pw_s[groups[-1][0]] + max(PW_ABS, PW_REL * pw_s[groups[-1][0]]):
            groups[-1].append(i)
        else:
            groups.append([i])
    emitters, left = [], 0
    for g in sorted(groups, key=lambda g: -len(g)):
        g = sorted(g, key=lambda i: toa_s[i])
        em, rest = _emitters([toa_s[i] for i in g], [pw_s[i] for i in g])
        for e in em:
            if e["n_pulses"] >= MIN_REPORT:
                emitters.append(e)
            else:
                left += e["n_pulses"]
        left += len(rest)
    emitters.sort(key=lambda e: -e["n_pulses"])
    return {"emitters": emitters, "unassigned": left}
