#!/usr/bin/env python3
"""
spectrum_survey.py -- what is on the air in a wideband SigMF capture, and
decode whatever has a decoder.

One pass over the capture, several independent detectors, one inventory:

  spurs     : CW tones / LO leakage from a median PSD (narrow, persistent)
  Wi-Fi     : wifi_detect (L-STF/L-LTF/L-SIG receiver, blind 1 MHz scan) +
              wifi_payload (non-HT MAC frames, FCS-verified) + undecodable
              Wi-Fi preambles
  15 kHz OFDM (LTE-like numerology: 66.67 us symbol + 4.7 us CP)
            : cyclic-prefix detector folded over the 71.4 us symbol period,
              run on an FFT-domain channelizer at 1.4 / 10 / 20 MHz widths.
              Wi-Fi (3.2 us symbols) does not correlate at 66.67 us, so this
              finds DJI OcuSync-class links (RC uplink hops, video downlink,
              DroneID) underneath busy Wi-Fi.
  DroneID   : 10 MHz / ~0.64 ms 15 kHz-OFDM bursts are handed to the DroneID
              decoder (droneid-decoder-reference/src/droneid_api.py) if found
  residual  : strong time-frequency energy no detector explains

Labels for non-Wi-Fi OFDM are inferred from numerology and burst shape, not
decoded; they say "-like" for that reason.

Usage
  python spectrum_survey.py cap.sigmf-meta [--out prefix] [--droneid-src DIR]
  -> prefix_inventory.json (emitters + detections), prefix_detections.csv,
     prefix.sigmf-meta (annotations), prefix.png (overlay)
"""
from __future__ import annotations
import json, os, sys
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict, field
import numpy as np
from scipy import ndimage

from spectrum_intel.core.spectral import Channelizer, find_spurs, _cp_metric, _folded_scores  # noqa: F401
from spectrum_intel.core.ofdm15k import (PROFILES, ofdm15k_scan, refine_ofdm15k, NUMEROLOGIES, ofdm_numerology, classify_numerology, label_ofdm15k, merge_ofdm15k, _cycle, track_hops, hop_summary)  # noqa: F401
from spectrum_intel.receivers.wifi import ofdm_rx as W
from spectrum_intel.receivers.wifi import payload as WP
# ============================================================== records
@dataclass
class Detection:
    kind: str                  # wifi_ppdu | wifi_undecodable | ofdm15k | spur | unclassified
    label: str                 # human label, e.g. "wifi_VHT", "ocusync_uplink_like"
    start_us: float
    end_us: float
    f_lo_mhz: float
    f_hi_mhz: float
    center_mhz: float
    snr_db: float | None = None
    decoded: bool = False      # carries content we actually demodulated (L-SIG, MAC frame, DroneID)
    info: dict = field(default_factory=dict)


# ============================================================== DroneID hook
def try_droneid(x, fs, fc, det, src_dir):
    if not src_dir or not os.path.isdir(src_dir):
        return dict(status="decoder_not_found")
    sys.path.insert(0, src_dir)
    try:
        from droneid_api import decode_ddc_packet
    except Exception as e:                                           # pragma: no cover
        return dict(status=f"import_failed: {e}")
    s = int(det.start_us * 1e-6 * fs); e = int(det.end_us * 1e-6 * fs)
    pad = int(50e-6 * fs)
    s, e = max(0, s - pad), min(len(x), e + pad)
    n = np.arange(s, e)
    iq = (x[s:e] * np.exp(-2j * np.pi * (det.center_mhz * 1e6 - fc) / fs * n)).astype(np.complex64)
    r = decode_ddc_packet(iq, sample_rate=fs, strip_pre_buffer=False)
    return dict(status=r.get("status"), success=r.get("success"),
                payload=(r.get("packet") or {}).get("payload"), error=r.get("error"))


# ============================================================== residual energy
def residual_energy(x, fs, fc, known, nfft=512, hop=256, thr_db=12.0, min_us=40, min_mhz=0.5,
                    pad_us=20, pad_mhz=1.5, top=50):
    """Strong time-frequency regions no detector accounts for."""
    n = (len(x) - nfft) // hop
    idx = np.arange(nfft)[None, :] + hop * np.arange(n)[:, None]
    P = np.abs(np.fft.fftshift(np.fft.fft(x[idx] * np.hanning(nfft), axis=1), axes=1)) ** 2
    Pd = 10 * np.log10(ndimage.uniform_filter(P, (3, 3)) + 1e-20)
    floor = ndimage.median_filter(np.percentile(Pd, 10, axis=0), 31)
    t_us = (np.arange(n) * hop + nfft / 2) / fs * 1e6
    f_mhz = ((np.arange(nfft) - nfft // 2) / nfft * fs + fc) / 1e6
    m = (Pd - floor) > thr_db
    for d in known:                                                   # blank everything explained
        ti = (t_us >= d.start_us - pad_us) & (t_us <= d.end_us + pad_us)
        fi = (f_mhz >= d.f_lo_mhz - pad_mhz) & (f_mhz <= d.f_hi_mhz + pad_mhz)
        m[np.ix_(ti, fi)] = False
    m = ndimage.binary_opening(m, np.ones((3, 3)))
    lab, k = ndimage.label(m)
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab)):
        t0, t1 = t_us[sl[0].start], t_us[sl[0].stop - 1]
        f0, f1 = f_mhz[sl[1].start], f_mhz[sl[1].stop - 1]
        if t1 - t0 < min_us or f1 - f0 < min_mhz:
            continue
        cells = lab[sl] == i + 1
        snr = float(np.percentile((Pd - floor)[sl][cells], 90))
        out.append(Detection("unclassified", "unclassified_energy", round(t0, 1), round(t1, 1),
                             round(f0, 3), round(f1, 3), round((f0 + f1) / 2, 3), round(snr, 1),
                             info=dict(fill=round(float(cells.mean()), 2))))
    out.sort(key=lambda d: -(d.end_us - d.start_us) * (d.f_hi_mhz - d.f_lo_mhz))
    return out[:top]


# ============================================================== survey
def survey(meta_path, data_path=None, out="survey", droneid_src=None, wifi_payloads=True, log=print):
    x, fs, fc, meta = W.load_sigmf(meta_path, data_path)
    dur_us = len(x) / fs * 1e6
    dets: list[Detection] = []
    log(f"{len(x)} samples @ {fs/1e6:g} MSps, fc {fc/1e6:g} MHz, {dur_us/1e3:.2f} ms")

    # ---- spurs
    spurs = find_spurs(x, fs, fc)
    for s in spurs:
        dets.append(Detection("spur", "lo_leakage" if "DC" in s["note"] else "cw_tone", 0.0, round(dur_us, 1),
                              s["freq_mhz"] - 0.02, s["freq_mhz"] + 0.02, s["freq_mhz"],
                              s["prominence_db"], info=s))
    log(f"spurs: {[(s['freq_mhz'], s['note']) for s in spurs]}")

    # ---- Wi-Fi
    det = W.WifiDetector(fs, fc, scan_step_hz=1e6)
    pkts = det.run(x)
    bursts = W.undecoded_wifi(det, x, pkts)
    for p in pkts:
        dets.append(Detection("wifi_ppdu", f"wifi_{p.fmt}", p.start_us, p.end_us, p.f_lo_mhz, p.f_hi_mhz,
                              p.center_mhz, p.snr_db, decoded=True,
                              info=dict(channel=p.channel, bw_mhz=p.bandwidth_mhz, rate=p.lsig_rate_mbps,
                                        length=p.lsig_length, mcs=p.mcs, guess=p.guess,
                                        off_grid=p.extra.get("off_grid"))))
    for b in bursts:
        dets.append(Detection("wifi_undecodable", b.label, b.start_us, b.start_us + b.duration_us,
                              b.f_lo_mhz, b.f_hi_mhz, b.center_mhz, None, info=dict(reason=b.reason)))
    frames = []
    if wifi_payloads:
        res = WP.run(meta_path, data_path, out + "_wifi_frames", scan=True)
        frames = [r for r in res if r.get("fcs_ok")]
    log(f"Wi-Fi: {len(pkts)} PPDUs {dict(Counter(p.fmt for p in pkts))}, "
        f"{len(bursts)} undecodable, {len(frames)} FCS-verified MAC frames")

    # ---- 15 kHz OFDM (OcuSync-class)
    notch = [(s["offset_hz"], 30e3) for s in spurs]
    ch = Channelizer(x.astype(np.complex128), fs, notch)
    hits, stats = ofdm15k_scan(ch, fc, usable_half=min(fs, 56e6) / 2)
    ob = merge_ofdm15k([b for b in (refine_ofdm15k(ch, h) for h in hits) if b is not None])
    ob = merge_ofdm15k(ob + track_hops(ch, fc, ob, stats, len(x) / fs, usable_half=min(fs, 56e6) / 2))
    o_dets = []
    for b in ob:
        lab, why = label_ofdm15k(b)
        cen = (fc + b["f_off"]) / 1e6
        d = Detection("ofdm15k", lab, round(b["start_s"] * 1e6, 1), round(b["end_s"] * 1e6, 1),
                      round(cen - b["bw_hz"] / 2e6, 3), round(cen + b["bw_hz"] / 2e6, 3), round(cen, 3),
                      None, info=dict(why=why, profile=b["profile"], cp_score=round(b["score"], 3),
                                      n_symbols=b["nsym"], bw_mhz=round(b["bw_hz"] / 1e6, 2),
                                      symbol_us=b["period_us"], cp_peaks=b["cp_peaks"],
                                      cp_scatter_us=b["cp_scatter_us"], tracked=b.get("tracked", False)))
        if lab == "droneid_candidate":
            d.info["droneid"] = try_droneid(x, fs, fc, d, droneid_src)
            d.decoded = bool(d.info["droneid"].get("success"))
        o_dets.append(d)
    dets += o_dets
    log(f"15 kHz OFDM bursts: {dict(Counter(d.label for d in o_dets))}")

    # ---- residual, labelled by blind numerology where an OFDM symbol length shows
    resid = residual_energy(x, fs, fc, dets)
    for d in resid:
        nm = ofdm_numerology(ch, d.center_mhz * 1e6 - fc, max((d.f_hi_mhz - d.f_lo_mhz) / 2 * 1e6, 0.5e6),
                             d.start_us * 1e-6, d.end_us * 1e-6)
        if nm:
            d.label, why = classify_numerology(nm[0])
            d.info.update(tu_us=round(nm[0] * 1e6, 2), cp_corr=round(nm[1], 3), why=why)
    dets += resid
    log(f"unclassified energy regions: {len(resid)}")

    dets.sort(key=lambda d: (d.start_us, d.center_mhz))
    emitters = build_emitters(dets, frames, fc, dur_us)
    inv = dict(capture=dict(file=os.path.basename(meta_path), fc_mhz=fc / 1e6, fs_msps=fs / 1e6,
                            duration_ms=round(dur_us / 1e3, 3),
                            description=meta["global"].get("core:description", ""),
                            datetime=meta["captures"][0].get("core:datetime")),
               emitters=emitters, detections=[asdict(d) for d in dets],
               wifi_frames=[{k: f.get(k) for k in ("start_us", "center_mhz", "subtype", "addr1_ra", "addr1",
                                                   "addr2", "addr2_ta", "bssid", "ssid", "seq",
                                                   "addr1_ra_vendor", "addr2_vendor", "bssid_vendor")
                             if f.get(k) is not None} for f in frames])
    return inv, x, fs, fc, meta, dets


def build_emitters(dets, frames, fc, dur_us):
    """Group detections into emitters/links a person can reason about."""
    em = []
    # Wi-Fi: one entry per occupied centre
    by_c = defaultdict(list)
    for d in dets:
        if d.kind == "wifi_ppdu":
            by_c[round(d.center_mhz)].append(d)
    for c, ds in sorted(by_c.items()):
        air = sum(d.end_us - d.start_us for d in ds)
        fr = [f for f in frames if abs(f["center_mhz"] - c) < 3]
        devs = {}
        for f in fr:
            for k in ("addr1_ra", "addr1", "addr2", "addr2_ta", "bssid"):
                m = f.get(k)
                if m and not (int(m[:2], 16) & 1):
                    devs[m] = f.get(k + "_vendor") or WP.vendor(m)       # control frames carry no vendor key
        em.append(dict(type="Wi-Fi", center_mhz=c, ppdus=len(ds),
                       formats=dict(Counter(d.label.replace("wifi_", "") for d in ds)),
                       bandwidths_mhz=sorted({d.info["bw_mhz"] for d in ds}),
                       off_grid=any(d.info.get("off_grid") for d in ds),
                       airtime_pct=round(100 * air / dur_us, 1),
                       devices=devs, ssids=sorted({f["ssid"] for f in fr if f.get("ssid")})))
    und = [d for d in dets if d.kind == "wifi_undecodable"]
    if und:
        em.append(dict(type="Wi-Fi (undecodable)", count=len(und),
                       centres_mhz=dict(Counter(round(d.center_mhz) for d in und)),
                       reasons=dict(Counter(d.info.get("reason", "").split(":")[0] for d in und))))
    # 15 kHz OFDM
    for lab in ("ocusync_uplink_like", "ocusync_video_like", "droneid_candidate", "ofdm15k_other"):
        ds = [d for d in dets if d.label == lab]
        if not ds:
            continue
        e = dict(type=lab, count=len(ds), why=ds[0].info["why"])
        if lab == "ocusync_uplink_like":
            e.update(hop_summary(ds, fc))
        if lab == "droneid_candidate":
            e["decoded"] = [d.info.get("droneid") for d in ds]
        em.append(e)
    for d in dets:
        if d.kind == "spur":
            em.append(dict(type=d.label, center_mhz=d.center_mhz, prominence_db=d.snr_db))
    res = [d for d in dets if d.kind == "unclassified"]
    if res:
        em.append(dict(type="unclassified energy", count=len(res), by_numerology=dict(Counter(d.label for d in res)),
                       largest=[dict(t_us=(d.start_us, d.end_us), f_mhz=(d.f_lo_mhz, d.f_hi_mhz), snr_db=d.snr_db,
                                     label=d.label) for d in res[:5]]))
    return em


# ============================================================== outputs
COLORS = {"wifi_ppdu": "#ff3b3b", "wifi_undecodable": "#ffa500", "ofdm15k": "#00e5ff",
          "spur": "#ffffff", "unclassified": "#b0b0b0"}


def plot(x, fs, fc, dets, path, nfft=1024):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    n = len(x) // nfft * nfft
    S = np.fft.fftshift(np.fft.fft(x[:n].reshape(-1, nfft) * np.hanning(nfft), axis=1), axes=1)
    S = 10 * np.log10(np.abs(S) ** 2 + 1e-12).T
    fig, ax = plt.subplots(figsize=(18, 8), dpi=110)
    ax.imshow(S, aspect="auto", origin="lower", cmap="viridis",
              extent=[0, n / fs * 1e3, (fc - fs / 2) / 1e6, (fc + fs / 2) / 1e6],
              vmin=np.percentile(S, 5), vmax=np.percentile(S, 99.7))
    for d in dets:
        if d.kind == "spur":
            continue
        ax.add_patch(Rectangle((d.start_us / 1e3, d.f_lo_mhz), (d.end_us - d.start_us) / 1e3,
                               d.f_hi_mhz - d.f_lo_mhz, fill=False,
                               lw=1.6 if d.kind == "ofdm15k" else 0.9, ec=COLORS[d.kind],
                               ls="--" if d.kind in ("wifi_undecodable", "unclassified") else "-"))
        if d.kind == "ofdm15k":
            ax.text(d.start_us / 1e3, d.f_hi_mhz + 0.3, d.label.replace("_like", ""), color=COLORS[d.kind], fontsize=7)
    for k, c in COLORS.items():
        if k != "spur":
            ax.plot([], [], color=c, label=k)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.6)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("frequency (MHz)")
    cnt = Counter(d.kind for d in dets)
    ax.set_title("Spectrum survey  |  " + "  ".join(f"{k}: {v}" for k, v in cnt.items()))
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def write_outputs(inv, x, fs, fc, meta, dets, out, make_plot=True):
    def clean(o):
        if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [clean(v) for v in o]
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.bool_): return bool(o)
        return o
    json.dump(clean(inv), open(out + "_inventory.json", "w"), indent=1)
    import csv
    with open(out + "_detections.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "label", "start_us", "end_us", "f_lo_mhz", "f_hi_mhz", "center_mhz", "snr_db",
                    "decoded", "info"])
        for d in dets:
            w.writerow([d.kind, d.label, d.start_us, d.end_us, d.f_lo_mhz, d.f_hi_mhz, d.center_mhz, d.snr_db,
                        d.decoded, json.dumps(clean(d.info))])
    m = json.loads(json.dumps(meta))
    m["annotations"] = sorted([{
        "core:sample_start": int(d.start_us * 1e-6 * fs),
        "core:sample_count": max(1, int((d.end_us - d.start_us) * 1e-6 * fs)),
        "core:freq_lower_edge": d.f_lo_mhz * 1e6, "core:freq_upper_edge": d.f_hi_mhz * 1e6,
        "core:label": d.label, "core:comment": json.dumps(clean(d.info))[:500],
    } for d in dets], key=lambda a: a["core:sample_start"])
    json.dump(m, open(out + ".sigmf-meta", "w"), indent=2)
    if make_plot:
        plot(x, fs, fc, dets, out + ".png")


def print_report(inv):
    c = inv["capture"]
    print(f"\n=== {c['file']}  fc {c['fc_mhz']} MHz  {c['fs_msps']} MSps  {c['duration_ms']} ms")
    if c["description"]:
        print(f"    {c['description']}")
    for e in inv["emitters"]:
        t = e["type"]
        if t == "Wi-Fi":
            devs = ", ".join(f"{m} ({v})" for m, v in e["devices"].items())
            print(f"  Wi-Fi  {e['center_mhz']} MHz{' (off-grid)' if e['off_grid'] else ''}: {e['ppdus']} PPDUs "
                  f"{e['formats']} BW {e['bandwidths_mhz']} MHz, airtime {e['airtime_pct']}%"
                  + (f"\n         devices: {devs}" if devs else "")
                  + (f"\n         SSIDs: {e['ssids']}" if e["ssids"] else ""))
        elif t == "Wi-Fi (undecodable)":
            print(f"  Wi-Fi preambles not decodable: {e['count']} at {e['centres_mhz']}")
        elif t == "ocusync_uplink_like":
            print(f"  OcuSync-like uplink (15 kHz OFDM, {e['bw_mhz']} MHz, {e['burst_us']:.0f} us bursts): "
                  f"{e['n_bursts']} hops at {e['centres_mhz']} MHz, intervals {e['intervals_ms']} ms"
                  + (f", cycle ~{e['cycle_ms']} ms" if "cycle_ms" in e else ""))
        elif t in ("ocusync_video_like", "droneid_candidate", "ofdm15k_other"):
            print(f"  {t}: {e['count']}  ({e['why']})")
            for dd in e.get("decoded", []):
                print(f"         DroneID: {dd}")
        elif t in ("lo_leakage", "cw_tone"):
            print(f"  {t} at {e['center_mhz']} MHz (+{e['prominence_db']} dB)")
        elif t == "unclassified energy":
            print(f"  unclassified energy regions: {e['count']}  by blind numerology: {e['by_numerology']}")
    if inv["wifi_frames"]:
        print(f"  FCS-verified Wi-Fi frames: {len(inv['wifi_frames'])}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Survey every signal in a SigMF capture; decode what we can")
    ap.add_argument("meta")
    ap.add_argument("--data")
    ap.add_argument("--out", default="survey")
    from spectrum_intel.receivers.dji._paths import droneid_src
    ap.add_argument("--droneid-src", default=droneid_src())
    ap.add_argument("--no-payload", action="store_true", help="skip Wi-Fi MAC frame decoding")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    inv, x, fs, fc, meta, dets = survey(a.meta, a.data, a.out, a.droneid_src, not a.no_payload)
    write_outputs(inv, x, fs, fc, meta, dets, a.out, not a.no_plot)
    print_report(inv)
