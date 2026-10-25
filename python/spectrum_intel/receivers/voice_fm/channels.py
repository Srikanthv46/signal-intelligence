"""Walkie-talkie channel plans and the sub-audible signalling tables (CTCSS tones, DCS codes)."""
from __future__ import annotations

# service -> list of (channel, nominal Hz); channel bandwidth Hz
PLANS: dict[str, tuple[list[tuple[int, float]], float]] = {
    "PMR446": ([(n, 446.00625e6 + 12.5e3 * (n - 1)) for n in range(1, 17)], 12.5e3),
    "LPD433": ([(n, 433.075e6 + 25e3 * (n - 1)) for n in range(1, 70)], 25e3),
    "FRS/GMRS": ([(n, f * 1e6) for n, f in zip(range(1, 23), (
        462.5625, 462.5875, 462.6125, 462.6375, 462.6625, 462.6875, 462.7125,
        467.5625, 467.5875, 467.6125, 467.6375, 467.6625, 467.6875, 467.7125,
        462.5500, 462.5750, 462.6000, 462.6250, 462.6500, 462.6750, 462.7000, 462.7250))], 12.5e3),
}
SPACING = {"PMR446": 12.5e3, "LPD433": 25e3, "FRS/GMRS": 12.5e3}

# EIA/TIA-603 CTCSS tones (Hz), plus the common 'non-EIA' extras radios offer
CTCSS = [67.0, 69.3, 71.9, 74.4, 77.0, 79.7, 82.5, 85.4, 88.5, 91.5, 94.8, 97.4, 100.0, 103.5, 107.2, 110.9,
         114.8, 118.8, 123.0, 127.3, 131.8, 136.5, 141.3, 146.2, 150.0, 151.4, 156.7, 159.8, 162.2, 165.5,
         167.9, 171.3, 173.8, 177.3, 179.9, 183.5, 186.2, 189.9, 192.8, 196.6, 199.5, 203.5, 206.5, 210.7,
         218.1, 225.7, 229.1, 233.6, 241.8, 250.3, 254.1]

# the 104 standard DCS codes (octal)
DCS_CODES = """023 025 026 031 032 036 043 047 051 053 054 065 071 072 073 074 114 115 116 122 125 131 132 134
143 145 152 155 156 162 165 172 174 205 212 223 225 226 243 244 245 246 251 252 255 261 263 265 266 271 274 306
311 315 325 331 332 343 346 351 356 364 365 371 411 412 413 423 431 432 445 446 452 454 455 462 464 465 466 503
506 516 523 526 532 546 565 606 612 624 627 631 632 654 662 664 703 712 723 731 732 734 743 754""".split()


def match_channel(f_hz: float, ppm: float = 0.0) -> dict | None:
    """Nearest walkie-talkie channel to a measured carrier. ``ppm`` is the receiver's
    frequency error (positive: its LO is high, so carriers read low); the tolerance is
    0.4 channel spacings or 12 ppm, whichever is smaller."""
    f = f_hz * (1 + ppm * 1e-6)
    best = None
    for svc, (chans, bw) in PLANS.items():
        tol = min(0.4 * SPACING[svc], 12e-6 * f)
        for n, nom in chans:
            if abs(f - nom) <= tol and (best is None or abs(f - nom) < abs(f - best[2])):
                best = (svc, n, nom, bw)
    if best is None:
        return None
    svc, n, nom, bw = best
    return {"service": svc, "channel": n, "nominal_mhz": round(nom / 1e6, 5), "channel_bw_khz": bw / 1e3,
            "offset_hz": round(f_hz - nom), "offset_ppm": round((f_hz - nom) / nom * 1e6, 2),
            "receiver_ppm_applied": ppm}
