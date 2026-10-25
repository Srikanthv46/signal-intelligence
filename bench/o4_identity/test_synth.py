#!/usr/bin/env python3
"""Bench: o4.py on synthetic O4-like links with known crystal offsets.

    cd bench/o4_identity && nice -n 19 ../../.venv/bin/python test_synth.py
"""
import numpy as np

import o4
from synth import o4_tx, pilot_table

NS = 1_000_000                       # one contiguous stretch (~16 ms)
FC = 5840e6
results = []


def check(name, ok, msg):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def link(eps, f_nom, snr, bw=20, seed=1, pilots=None, ramp=0.0, t0=12345.6, fc=FC):
    return o4_tx(NS, bw_mhz=bw, eps_ppm=eps, f_off_hz=(f_nom - fc) + f_nom * eps * 1e-6, snr_db=snr, seed=seed,
                 pilots=pilots, amp_ramp_db=ramp, t0=t0)


def measure(x, f_nom, bw=20, f0_err=30e3, fc=FC):
    L = o4.analyse_link(x, fc, f_nom + f0_err, bw * 0.9e6, 0, NS)   # detector-like centre error
    rc = o4.resolve_carrier(L["f0"], L["cfo"], L["sco"]["ppm"] or 0.0)
    car, ppm = o4.carrier_ppm(L["f0"], L["cfo"], rc["n"], rc["f_nom"])
    return L, rc, ppm


# 1. clock and carrier ppm vs SNR (20 MHz, 5823.5 MHz, eps = +1.5 ppm)
for snr in (20, 10, 5):
    x = link(1.5, 5823.5e6, snr)
    L, rc, ppm = measure(x, 5823.5e6)
    s = L["sco"]
    ok = s["ppm"] is not None and abs(s["ppm"] - 1.5) < max(0.1, 4 * s["sigma"]) and rc["f_nom"] == 5823.5e6 \
        and abs(ppm - 1.5) < 0.02
    check(f"ppm @ {snr} dB", ok, f"clock {s['ppm']:+.3f} +- {s['sigma']:.3f} ppm (n={s['n']}), nominal "
          f"{rc['f_nom'] / 1e6:.3f} MHz (raster residual {rc['raster_residual_hz']:+.0f} Hz), carrier {ppm:+.4f} ppm "
          f"[true +1.5], grid quality {L['grid_quality']:.1f}")

# 2. power ramp of 20 dB across the stretch (approach / recede): no effect
x = link(-2.3, 5823.5e6, 15, ramp=20.0)
L, rc, ppm = measure(x, 5823.5e6)
check("20 dB power ramp", abs(ppm + 2.3) < 0.02 and abs(L["sco"]["ppm"] + 2.3) < 0.15,
      f"carrier {ppm:+.4f} ppm, clock {L['sco']['ppm']:+.3f} [true -2.3]")

# 3. 10 MHz mode, same crystal
x10 = link(1.5, 5730.5e6, 10, bw=10, fc=5728e6)
L10, rc10, ppm10 = measure(x10, 5730.5e6, bw=10, fc=5728e6)
check("10 MHz mode ppm", abs(ppm10 - 1.5) < 0.03 and rc10["f_nom"] == 5730.5e6,
      f"carrier {ppm10:+.4f} ppm, clock {L10['sco']['ppm']:+.3f} +- {L10['sco']['sigma']:.3f}, nominal "
      f"{rc10['f_nom'] / 1e6:.3f} MHz")

# 4. sequence identity: 10 MHz vs 20 MHz of the same pilot table -> one common shift; another table -> none
x20 = link(1.5, 2415.5e6, 20, fc=2428e6)
L20 = o4.analyse_link(x20, 2428e6, 2415.5e6, 18e6, 0, NS)
L10 = o4.analyse_link(x10, 5728e6, 5730.5e6, 9e6, 0, NS)
b, s = o4.match(L10["rows"], L10["K"], L20["rows"], L20["K"], 200)
st = o4.match_stats(b, s, 200)
check("10 vs 20 MHz same sequence", st["p"] < 1e-6 and st["shift"] == -150,
      f"shift {st['shift']} x{st['hits']}/{st['n']} (p {st['p']:.2g}) [expected -150]")
xo = link(1.5, 2415.5e6, 20, fc=2428e6, pilots=pilot_table(seed=99), seed=5)
Lo = o4.analyse_link(xo, 2428e6, 2415.5e6, 18e6, 0, NS)
b, s = o4.match(Lo["rows"], Lo["K"], L20["rows"], L20["K"], 200)
st = o4.match_stats(b, s, 200)
check("different sequence -> no match", st["p"] > 1e-4, f"best shift {st['shift']} x{st['hits']}/{st['n']} (p {st['p']:.2g})")

# 5. two transmitters of one model, 0.2 ppm apart, on two channels of one capture: resolved
xa = link(1.5, 5823.5e6, 15, seed=11, t0=1000.3)
xb = link(1.7, 5845.5e6, 15, seed=12, t0=20000.7)
x2 = xa + xb
La, rca, pa = measure(x2, 5823.5e6)
Lb, rcb, pb = measure(x2, 5845.5e6)
check("two drones 0.2 ppm apart", abs((pb - pa) - 0.2) < 0.03,
      f"A {pa:+.4f}, B {pb:+.4f} ppm, difference {pb - pa:+.4f} [true +0.200]")

print(f"\n{sum(results)}/{len(results)} passed")
