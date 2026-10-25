"""Rule-based classifier: India band prior + measured features -> label.

Every rule returns (label, score, reasons) or None. Scores add up the
independent evidence a rule found (band match, bandwidth, numerology, timing),
so a label backed by more independent features scores higher. The best rule
wins; below MIN_CONF the signal is reported as "unknown" with the candidates
listed. Labels that are inferred (not decoded) say "-like".

A decoder result always overrides the rules (its label is backed by a CRC /
correlation against the standard's sequences).
"""
from __future__ import annotations

from ..model import Signal
from ..core import bands_india as B

MIN_CONF = 0.5
MHz, kHz = 1e6, 1e3


def _near(v, targets, tol):
    return v is not None and any(abs(v - t) <= tol * t for t in targets)


def _keys(s):
    return set(B.band_keys(s.f0, s.f1)) | {b[2] for b in B.lookup(s.fc)}


def r_cw(s, f):
    if s.modulation.startswith("CW"):
        return ("CW carrier (unmodulated tone)-like", 0.55,
                [f"bandwidth {s.bw:.0f} Hz: unmodulated carrier or receiver spur (sweep check decides)"])


def r_wifi_ofdm(s, f):
    k = _keys(s)
    if not k & {"ism_24", "unii_12", "unii_2e", "ism_58", "its_59", "wifi6e"}:
        return None
    sc, why = 0.25, ["inside a delicensed Wi-Fi band"]
    if _near(f.get("ofdm_scs_hz"), [312.5e3], 0.05):
        sc += 0.35; why.append("OFDM subcarrier spacing 312.5 kHz (802.11a/g/n/ac)")
    elif f.get("ofdm_scs_hz"):
        return None
    if _near(s.bw, [16.6e6, 17.5e6, 36e6, 38e6, 76e6, 78e6], 0.12):
        sc += 0.15; why.append(f"occupied bandwidth {s.bw / MHz:.1f} MHz matches 20/40/80 MHz channel")
    if not f.get("continuous") and f.get("n_bursts", 0) >= 1:
        sc += 0.05; why.append("bursty (CSMA)")
    return ("Wi-Fi 802.11 OFDM-like", sc, why) if sc >= 0.5 else None


def r_wifi_dsss(s, f):
    if "ism_24" in _keys(s) and s.modulation.startswith("DSSS"):
        return ("Wi-Fi 802.11b DSSS-like", 0.75, ["2.4 GHz ISM", "1 Msym/s Barker DSSS clock"])


def _cell_band(k):
    dl = {"b28_dl": "B28/n28", "b5_dl": "B5/n5", "b8_dl": "B8/n8", "b3_dl": "B3/n3", "b1_dl": "B1/n1"}
    tdd = {"b40": "B40/n40", "b41": "B41/n41", "n78": "n78"}
    ul = {"b28_ul": "B28", "b5_ul": "B5", "b8_ul": "B8", "b3_ul": "B3", "b1_ul": "B1"}
    for d, kind in ((dl, "dl"), (tdd, "tdd"), (ul, "ul")):
        for key, name in d.items():
            if key in k:
                return kind, name
    return None, None


def r_cellular_ofdm(s, f):
    kind, band = _cell_band(_keys(s))
    if not kind:
        return None
    scs = f.get("ofdm_scs_hz")
    sc, why = 0.25, [f"inside India IMT band {band} ({kind.upper()})"]
    tech = None
    if _near(scs, [15e3], 0.04):
        sc += 0.3; why.append("OFDM SCS 15 kHz (LTE / NR-15)"); tech = "LTE"
    elif _near(scs, [30e3], 0.04):
        sc += 0.3; why.append("OFDM SCS 30 kHz (NR)"); tech = "5G NR"
    elif scs:
        return None
    if _near(s.bw, [1.4e6, 3e6, 5e6, 10e6, 15e6, 20e6, 40e6, 50e6, 60e6, 80e6, 100e6], 0.12):
        sc += 0.1; why.append(f"bandwidth {s.bw / MHz:.1f} MHz on a 3GPP channel size")
    per = f.get("burst_period_s")
    if kind == "dl" and f.get("continuous"):
        sc += 0.1; why.append("always on (downlink reference signals)")
    if kind == "tdd" and per and any(abs(per - p) < 0.1e-3 for p in (2.5e-3, 5e-3, 10e-3)):
        sc += 0.15; why.append(f"TDD frame periodicity {per * 1e3:.2f} ms")
    if kind == "ul" and not f.get("continuous"):
        sc += 0.05; why.append("bursty (UE uplink scheduling)")
    if tech is None:
        tech = "LTE/NR"
        if sc < 0.45:
            return None
    lab = {"dl": f"{tech}-like downlink {band}", "tdd": f"{tech}-like TDD {band}",
           "ul": f"{tech}-like uplink {band}"}[kind]
    return lab, sc, why


def r_gsm(s, f):
    k = _keys(s)
    if not k & {"b8_dl", "b3_dl", "b8_ul", "b3_ul"}:
        return None
    if not 120e3 <= s.bw <= 320e3:
        return None
    sc, why = 0.35, ["GSM band", f"bandwidth {s.bw / kHz:.0f} kHz (GSM 200 kHz carrier)"]
    if "constant" in s.modulation or "GMSK" in s.modulation:
        sc += 0.2; why.append("constant envelope (GMSK)")
    if f.get("continuous") and k & {"b8_dl", "b3_dl"}:
        sc += 0.1; why.append("continuous (BCCH carrier)")
    arfcn = B.gsm_arfcn(s.fc)
    return (f"GSM-like {'downlink' if k & {'b8_dl', 'b3_dl'} else 'uplink'}" +
            (f" ARFCN {arfcn}" if arfcn is not None else ""), sc, why)


def r_umts(s, f):
    if "b1_dl" in _keys(s) and 3.5e6 <= s.bw <= 5.2e6 and not f.get("ofdm_scs_hz") and f.get("continuous"):
        sc, why = 0.5, ["B1 downlink", f"bandwidth {s.bw / MHz:.2f} MHz (3.84 Mcps WCDMA)", "continuous, no OFDM CP"]
        if _near(f.get("cyclic_line_hz"), [3.84e6], 0.01):
            sc += 0.2; why.append("3.84 MHz chip-rate cyclic line")
        return "UMTS/WCDMA-like downlink B1", sc, why


def r_ble(s, f):
    if "ism_24" not in _keys(s) or not 0.6e6 <= s.bw <= 2.6e6:
        return None
    if f.get("continuous") or f.get("burst_dur_s_median", 1) > 2.2e-3:
        return None
    sc, why = 0.35, ["2.4 GHz ISM", f"bandwidth {s.bw / MHz:.2f} MHz, bursts {f['burst_dur_s_median'] * 1e6:.0f} us"]
    if "FSK" in s.modulation or "constant" in s.modulation:
        sc += 0.2; why.append("GFSK / constant envelope")
    for ch, fc in B.BLE_ADV.items():
        if abs(s.fc - fc) < 0.3e6:
            sc += 0.1; why.append(f"centred on BLE advertising channel {ch}")
    return "Bluetooth (BLE/BR)-like", sc, why


def r_adsb(s, f):
    if abs(s.fc - 1090e6) < 1.5e6 and f.get("burst_dur_s_median", 1) < 150e-6:
        return ("Mode S / ADS-B (1090 MHz)-like", 0.7,
                ["1090 MHz SSR reply frequency", f"short bursts ({f['burst_dur_s_median'] * 1e6:.0f} us)"])
    if abs(s.fc - 1030e6) < 1.5e6 and f.get("burst_dur_s_median", 1) < 60e-6:
        return "SSR interrogation (1030 MHz)-like", 0.65, ["1030 MHz SSR interrogation", "short pulses"]


PULSE_TB_MIN = 1.0              # bandwidth x pulse width of a real pulse
DME_PAIR_S = (10e-6, 45e-6)     # pulse-pair group width (two 3.5 us pulses, 12 us X / 30 us Y mode)


def r_pulsed(s, f):
    k = _keys(s)
    if not k & B.PULSED_BANDS:
        return None
    if any(abs(s.fc - f0) < 1.5e6 for f0 in (1030e6, 1090e6)):
        return None                    # SSR interrogations / Mode S replies (r_adsb), not DME
    # pulse widths from the 1 us envelope when measured (spectrogram burst boxes are far longer than the pulses)
    pulses = (f.get("pulse_n") or 0) >= 2
    dur = f["pulse_width_s_median"] if pulses else f.get("burst_dur_s_median", 1)
    if f.get("continuous") or dur > 200e-6:
        return None
    if (f.get("pulse_n") if pulses else f.get("n_bursts") or 0) < 2:
        return None                    # DME replies come in pulse pairs, radars repeat on their PRI: one blip is not
    if pulses and s.bw * dur < PULSE_TB_MIN:
        return None                    # a pulse of width tau occupies >~ 1/tau: shorter "pulses" are noise peaks
    if "unii_2e" in k and not (f.get("pri_s") or f.get("burst_period_s")):
        return None                    # shared with Wi-Fi: bursts there are packets unless they repeat on a PRI (DFS radar)
    dme = "arns_l" in k or "rnss_l5" in k
    if dme and pulses and not DME_PAIR_S[0] <= dur <= DME_PAIR_S[1]:
        return None                    # a DME reply is a pulse pair, 12 (X) or 30 (Y) us apart: one group 15-35 us
    name, radar_class, extra = _radar_class(s, f, k, dme, pulses, dur)
    sc = 0.5
    why = [B.describe(s.fc), f"short pulses ({dur * 1e6:.1f} us" +
           (f", {f['pulse_n']} pulse groups, 1 us envelope)" if pulses else ")")] + extra
    if f.get("pri_s"):                 # intervals within +-15 % of the shortest cluster
        sc += 0.15; why.append(f"PRI {f['pri_s'] * 1e6:.0f} us" +
                               (f" (staggered, spread {f['pri_spread'] * 100:.0f} %)" if f["pri_spread"] > 0.02 else ""))
    elif f.get("burst_period_s"):
        sc += 0.15; why.append(f"regular PRI {f['burst_period_s'] * 1e6:.0f} us")
    if f.get("lfm_sweep_hz"):
        sc += 0.1; why.append(f"LFM chirp {f['lfm_sweep_hz'] / 1e3:+.0f} kHz over {f['lfm_pulse_s'] * 1e6:.0f} us "
                              "(pulse compression)")
    if f.get("pulse_emitters"):        # core.deinterleave (research, signal gaps S10): evidence only
        why.append(f"deinterleaved into {len(f['pulse_emitters'])} emitter(s): "
                   + "; ".join(_emitter_text(e) for e in f["pulse_emitters"]))
    if radar_class:
        sc = min(0.8, sc + 0.1)
    return name, sc, why


def _emitter_text(e):
    st = e["stagger_us"]
    return (f"PRI {e['pri_us']:.1f} us" + (f" staggered {'/'.join(f'{v:.0f}' for v in st)} us" if st else "")
            + (f", jitter {e['jitter'] * 100:.1f} %" if e["jitter"] >= 0.01 else "")
            + f", {e['n_pulses']} pulses, PW {e['pw_us']:.1f} us")


def _radar_class(s, f, k, dme, pulses, dur):
    """Radar class from the India band plan, pulse width, PRI and LFM (research, signal gaps S10; file 05 sec. 4.1).
    Returns (label, radar_class, extra evidence); the class only raises the confidence, the label and evidence carry
    it. Every label keeps the '-like' suffix."""
    if dme:
        return "DME/TACAN pulse pairs-like", None, []
    pri = f.get("pri_s")
    fc = s.fc
    lfm = bool(f.get("lfm_sweep_hz"))
    long_lfm = lfm and (f.get("lfm_pulse_s") or 0) >= 10e-6
    if "radar_s_atc" in k:             # 2700-2900 MHz: ATC terminal radar and S-band weather radar share the band
        if long_lfm:                   # a long LFM (pulse-compression) pulse is the solid-state ATC signature
            return ("Pulsed radar (S-band ATC, solid-state)-like", "atc_s_solid_state",
                    ["S-band ATC band (2.7-2.9 GHz), long LFM pulse (pulse compression)"])
        if pulses and dur > 1.5e-6 and pri and 0.77e-3 <= pri <= 4e-3:
            return "Pulsed radar (weather, S-band)-like", "weather_s", ["S-band weather-radar band (2.7-2.9 GHz)"]
        return "Pulsed radar (S-band ATC)-like", "atc_s", ["S-band ATC band (2.7-2.9 GHz)"]
    if "radar_s_marine" in k:          # 2900-3100 MHz
        if dur <= 1.2e-6:
            return "Pulsed radar (marine S-band)-like", "marine_s", ["maritime navigation-radar band (2.9-3.1 GHz)"]
        return "Pulsed radar (S-band)-like", "s_band", ["S-band radar (2.9-3.1 GHz)"]
    if 5570e6 <= fc <= 5660e6:         # C-band weather radar (IMD DWR, 5600-5650 MHz)
        return "Pulsed radar (weather, C-band)-like", "weather_c", ["C-band weather-radar band (5.6 GHz)"]
    if "radar_5400" in k or "unii_2e" in k:
        return "Pulsed radar (C-band)-like", "c_band", ["C-band radiolocation (5.25-5.85 GHz)"]
    if ("rnss_l2_radar" in k or "arns_radar" in k or "fixed_radar_1350" in k):
        if 1270e6 <= fc <= 1295e6 and pri and pri <= 250e-6:
            return "Wind profiler radar-like", "wind_profiler", ["L-band wind-profiler range (~1280 MHz)"]
        if pulses and dur >= 1e-6 and pri and 1e-3 <= pri <= 4e-3:
            return ("Pulsed radar (L-band surveillance)-like", "surveillance_l",
                    ["L-band surveillance-radar band (1.21-1.4 GHz)"])
        return "Pulsed radar (L-band)-like", "l_band", ["L-band radar (1.21-1.4 GHz)"]
    if 440e6 <= fc <= 450e6 and pri and pri <= 250e-6:
        return "Wind profiler radar-like", "wind_profiler", ["UHF wind-profiler range (440-450 MHz)"]
    return "Pulsed radar-like", None, []


def r_dvbt2(s, f):
    if "uhf_tv" not in _keys(s) and "uhf_bcast_imt" not in _keys(s):
        return None
    if 6.5e6 <= s.bw <= 8e6 and f.get("continuous"):
        sc, why = 0.45, ["UHF TV band", f"bandwidth {s.bw / MHz:.2f} MHz (8 MHz raster)", "continuous"]
        if _near(f.get("ofdm_scs_hz"), [279e3 / 64 * 1, 1116, 2232, 4464, 558, 279], 0.05):
            sc += 0.35; why.append(f"OFDM SCS {f['ofdm_scs_hz']:.0f} Hz (DVB-T2 FFT mode)")
        elif f.get("ofdm_scs_hz"):
            sc += 0.15; why.append(f"OFDM SCS {f['ofdm_scs_hz']:.0f} Hz")
        return "DVB-T2 TV (Doordarshan)-like", sc, why


def r_narrowband(s, f):
    k = _keys(s)
    if not k & {"govt_mobile_300", "govt_mobile_335", "tetra_ppdr", "govt_mobile_406", "pmr_450",
                "govt_mobile_435", "trunk_800", "trunk_850", "metaids"}:
        return None
    if not 3e3 <= s.bw <= 30e3:
        return None
    if "tetra_ppdr" in k and 18e3 <= s.bw <= 27e3:
        return "TETRA-like (25 kHz pi/4-DQPSK)", 0.55, ["380-400 MHz public-safety band", "~25 kHz"]
    if "metaids" in k:
        return "Radiosonde-like telemetry", 0.5, ["400-406 MHz MetAids", f"{s.bw / kHz:.1f} kHz"]
    m = "analog FM / FSK voice-data" if "constant" in s.modulation or "FSK" in s.modulation else "narrowband"
    return f"PMR land-mobile ({m})-like", 0.55, [B.describe(s.fc), f"{s.bw / kHz:.1f} kHz channel"]


def r_srd(s, f):
    k = _keys(s)
    if k & {"ism_433", "ism_865"} and s.bw < 600e3 and not f.get("continuous"):
        return ("ISM short-range device (433/865 MHz)-like", 0.55,
                [B.describe(s.fc), f"bursty, {s.bw / kHz:.0f} kHz"])


def r_iridium(s, f):
    if "mss_iridium" in _keys(s) and 25e3 <= s.bw <= 60e3 and not f.get("continuous"):
        # Iridium bursts are DE-QPSK / BPSK at 25 ksym/s: without a PSK signature (x^2 / x^4 carrier line, or the
        # 25 kHz symbol line) a narrow burst here is as likely a wandering clock-comb line (live 1620.48 MHz)
        psk = [w for w, ok in (("x^2 carrier line", f.get("x2_line_ratio")), ("x^4 carrier line", f.get("x4_line_ratio")),
                               ("25 kHz symbol line", _near(f.get("cyclic_line_hz"), [25e3], 0.05))) if ok]
        if not psk:
            return None
        return ("Iridium-like TDMA burst", 0.55,
                ["1616-1626.5 MHz Iridium band", f"{s.bw / kHz:.0f} kHz burst", " + ".join(psk) + " (PSK)"])


RULES = [r_cw, r_wifi_ofdm, r_wifi_dsss, r_cellular_ofdm, r_gsm, r_umts, r_ble, r_adsb, r_pulsed,
         r_dvbt2, r_narrowband, r_srd, r_iridium]


def classify(s: Signal) -> None:
    f = s.features
    if "lo_leakage" in s.flags:
        s.label, s.confidence, s.modulation = "LO leakage (receiver artifact)", 0.95, "CW"
        s.evidence = "narrow tone exactly at the tuned centre frequency (DC bin): AD936x LO feed-through"
        return
    if "receiver_spur" in s.flags:
        fz = f.get("spur_freq_hz", s.fc)
        if f.get("spur_kind") == "sample_clock_harmonic":
            k = f.get("spur_harmonic")
            s.label, s.confidence = f"Receiver spur (sample-clock harmonic {k} x {f.get('spur_sample_rate_hz', 0) / MHz:g} MHz)", 0.9
            s.evidence = (f"narrow line at {fz / MHz:.6f} MHz = {k} x the sample rate within "
                          f"{abs(fz - k * f.get('spur_sample_rate_hz', 0)):.0f} Hz: generated by the receiver's own clock")
        else:
            s.label, s.confidence = "Receiver spur (locked to the receiver reference)", 0.8
            s.evidence = (f"narrow line at {fz / MHz:.6f} MHz, on the integer-MHz grid of the receiver's own frequency "
                          f"scale although the receiver error ({f.get('spur_receiver_ppm')} ppm) would move an "
                          f"external line {abs(f.get('spur_receiver_ppm', 0)) * fz * 1e-6 / 1e3:.1f} kHz off it")
        s.modulation = "CW"
        return
    if "iq_image" in s.flags:
        src = f.get("image_of_f_center_hz")
        s.label, s.confidence = "IQ image (receiver artifact)", 0.8
        s.evidence = (f"mirror of a >=15 dB stronger signal at {src / MHz:.3f} MHz about the LO, same burst "
                      f"timing: AD936x IQ-imbalance image")
        return
    cands = []
    for r in RULES:
        c = r(s, f)
        if c:
            cands.append(c)
    cands.sort(key=lambda c: -c[1])
    s.features["classifier_candidates"] = [{"label": c[0], "score": round(c[1], 2)} for c in cands[:4]]
    band = B.describe(s.fc)
    if cands and cands[0][1] >= MIN_CONF:
        lab, sc, why = cands[0]
        s.label, s.confidence = lab, round(min(sc, 0.9), 2)
        s.evidence = "; ".join(why + s.classify_trace[:2])
    else:
        s.label, s.confidence = "unknown", round(cands[0][1], 2) if cands else 0.2
        s.evidence = (f"band prior: {band}. Features: {s.modulation}, bw {s.bw / kHz:.0f} kHz, "
                      f"duty {f.get('duty_cycle')}" +
                      (f"; best candidate '{cands[0][0]}' only {cands[0][1]:.2f}" if cands else
                       "; no rule matched") + ("; " + "; ".join(s.classify_trace[:2]) if s.classify_trace else ""))
