"""Which transmitter sent each 802.11 PPDU -- including the ones whose MAC header is out of reach.

    attribute(recs, sifs_us=16.0)  -> transmitters (list of dicts); sets rec["ta"], rec["ta_evidence"]

`recs` are the pipeline's PPDU records, in any order: t0/t1 (s), fmt, `mac` (dissected MPDU when the
FCS passed), scrambler_seed / seed_ok, power_db, cfo_khz. Evidence, strongest first:

  1. "mac header": an FCS-valid MPDU names its transmitter (addr2 / TA).
  2. "response": a control response follows its frame after SIFS and names, in its RA, the station it
     answers. So the frame before an ACK / Block Ack / CTS addressed to X was sent by X, and an
     ACK / CTS (which has no TA) was sent by whoever the previous frame addressed.
  3. "scrambler counter": many chipsets (Atheros/Qualcomm among them) start each PPDU's scrambler
     from a counter that steps by one per transmission, 1..127. A transmitter whose known seeds step
     like that is "sequential". Between two of its PPDUs with seeds s_a and s_b, exactly
     (s_b - s_a) mod 127 transmissions happened; if that equals the number of PPDUs in between that are
     not someone else's, all of them are its (this attributes PPDUs no decoder can read, such as two-
     stream MIMO frames). A single PPDU whose own seed fits the counter on both sides is its too
     (allowing up to two transmissions the capture missed); at the ends of the counter only the next
     one or two steps count.
     Seeds are read from the SERVICE field, which decodes (and self-checks: 9 bits must descramble to
     zero) long before the payload does.
  4. Received power and carrier offset are reported per transmitter, and for a PPDU still unattributed
     the one transmitter within 3 dB is given as `ta_likely` (weak: power says where, not who).
"""
from __future__ import annotations

from collections import Counter

import numpy as np

RESPONSES = {"ack", "cts", "block-ack"}


def _nxt(seed, n=1):
    return (seed - 1 + n) % 127 + 1


def _delta(a, b):
    """Counter steps from seed a to seed b (1..127)."""
    return (b - a) % 127 or 127


def _set(r, ta, how):
    if r.get("ta") is None and ta:
        r["ta"], r["ta_evidence"] = ta, how
        return True
    return False


def attribute(recs, sifs_us=16.0):
    """Attribute every PPDU it can; return one summary row per transmitter."""
    recs = sorted(recs, key=lambda r: r["t0"])
    for r in recs:
        r.pop("ta", None), r.pop("ta_evidence", None), r.pop("ta_likely", None)
        m = r.get("mac") or {}
        r["_sub"] = m.get("subtype")
        r["_ra"] = m.get("ra")
        if m.get("ta"):
            _set(r, m["ta"], "mac header")
    # ---- 2. control responses (gap from the end of one PPDU to the start of the next)
    # sifs_us is the link's SIFS (16 us at a 20 MHz clock, scaled with the clock). Vendors stretch the
    # response delay for range (ACK timeout / coverage class): 11 MHz Microhard links answer 46-61 us
    # after the frame, so the window is generous; nothing else may start inside it (NAV).
    lo, hi = -4e-6, 84e-6 * sifs_us / 16.0
    for a, b in zip(recs, recs[1:]):
        gap = b["t0"] - a["t1"]
        if not lo <= gap <= hi:
            continue
        if b["_sub"] in RESPONSES and b["_ra"] and a["_sub"] not in RESPONSES:
            _set(a, b["_ra"], "response")
        if b["_sub"] in ("ack", "cts") and a["_ra"]:
            _set(b, a["_ra"], "response")
    # ---- 3. scrambler counters, until nothing changes
    sequential = {}
    for _ in range(4):
        changed = False
        by = {}
        for i, r in enumerate(recs):
            if r.get("ta") and r.get("seed_ok") and r.get("scrambler_seed"):
                by.setdefault(r["ta"], []).append(i)
        for ta, idx in by.items():
            pairs = list(zip(idx, idx[1:]))

            def between(i, j):
                return [k for k in range(i + 1, j) if recs[k].get("ta") in (None, ta)]

            good = [_delta(recs[i]["scrambler_seed"], recs[j]["scrambler_seed"]) <= len(between(i, j)) + 1
                    for i, j in pairs]
            sequential[ta] = len(pairs) >= 2 and np.mean(good) >= 0.8
            if not sequential[ta]:
                continue
            for i, j in pairs:                                      # exact gap fill
                bt = between(i, j)
                if bt and _delta(recs[i]["scrambler_seed"], recs[j]["scrambler_seed"]) == len(bt) + 1:
                    for k in bt:
                        changed |= _set(recs[k], ta, "scrambler counter")
        # single PPDUs whose own seed fits one sequential counter on both sides
        for k, r in enumerate(recs):
            if r.get("ta") or not (r.get("seed_ok") and r.get("scrambler_seed")):
                continue
            fits = []
            for ta, idx in by.items():
                if not sequential.get(ta):
                    continue
                before = [i for i in idx if i < k]
                after = [j for j in idx if j > k]
                both = bool(before and after)
                # Two-sided: allow up to two transmissions the capture did not resolve (a collision, a
                # clipped burst) -- a random seed fits both windows with probability < (5/127)^2.
                # One-sided: only the next step or two of a long counter.
                ok = bool(before or after)
                if before:
                    i = before[-1]
                    n = sum(1 for q in range(i + 1, k) if recs[q].get("ta") in (None, ta))
                    ok &= _delta(recs[i]["scrambler_seed"], r["scrambler_seed"]) <= (min(n + 3, 5) if both
                                                                                     else min(n + 1, 2))
                if after:
                    j = after[0]
                    n = sum(1 for q in range(k + 1, j) if recs[q].get("ta") in (None, ta))
                    ok &= _delta(r["scrambler_seed"], recs[j]["scrambler_seed"]) <= (min(n + 3, 5) if both
                                                                                     else min(n + 1, 2))
                ok &= both or len(idx) >= 4                     # one-sided fits need a long counter
                if ok:
                    fits.append(ta)
            if len(fits) == 1:
                changed |= _set(r, fits[0], "scrambler counter")
        if not changed:
            break
    # ---- 4. per-transmitter physical summary, and a weak power hint for what is left
    tx = {}
    for r in recs:
        if r.get("ta"):
            tx.setdefault(r["ta"], []).append(r)
    rows = []
    for ta, rs in sorted(tx.items(), key=lambda kv: -len(kv[1])):
        pw = [r["power_db"] for r in rs if r.get("power_db") is not None]
        cfo = [r["cfo_khz"] for r in rs if r.get("cfo_khz") is not None]
        seeds = [r["scrambler_seed"] for r in rs if r.get("seed_ok") and r.get("scrambler_seed")]
        vendor = next((r["mac"].get("ta_vendor") or r["mac"].get("addr2_ta_vendor") or r["mac"].get("addr2_vendor")
                       for r in rs if r.get("mac") and r["mac"].get("ta") == ta), None)
        rows.append({"ta": ta, "vendor": vendor, "n_ppdus": len(rs),
                     "evidence": dict(Counter(r["ta_evidence"] for r in rs)),
                     "formats": dict(Counter(r["fmt"] for r in rs)),
                     "frame_types": dict(Counter(r["_sub"] for r in rs if r["_sub"])),
                     "power_db_median": round(float(np.median(pw)), 1) if pw else None,
                     "power_db_spread": round(float(np.ptp(pw)), 1) if len(pw) > 1 else None,
                     "cfo_khz_median": round(float(np.median(cfo)), 2) if cfo else None,
                     "scrambler_counter": ("sequential" if sequential.get(ta) else
                                           "not sequential" if len(seeds) >= 3 else None),
                     "seeds": seeds[:64],
                     "t_first_ms": round(rs[0]["t0"] * 1e3, 3), "t_last_ms": round(rs[-1]["t0"] * 1e3, 3)})
    levels = {row["ta"]: row["power_db_median"] for row in rows if row["power_db_median"] is not None}
    for r in recs:
        if r.get("ta") or r.get("power_db") is None or not levels:
            continue
        near = [ta for ta, p in levels.items() if abs(p - r["power_db"]) <= 3]
        far = all(abs(p - r["power_db"]) > 6 for ta, p in levels.items() if ta not in near)
        if len(near) == 1 and far:
            r["ta_likely"] = near[0]
    for r in recs:
        r.pop("_sub", None), r.pop("_ra", None)
    unattributed = [r for r in recs if not r.get("ta")]
    for row in rows:
        row["n_likely_by_power"] = sum(1 for r in unattributed if r.get("ta_likely") == row["ta"])
    return rows
