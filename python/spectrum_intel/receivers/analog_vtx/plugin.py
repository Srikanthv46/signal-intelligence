"""Pipeline plugin: analog FM video transmitters (PAL/NTSC), demodulated to a picture (wraps .detector)."""
from __future__ import annotations


from spectrum_intel.core.bands_india import ANALOG_VTX_RANGES
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean, _in_ranges
from spectrum_intel.decode.base import result


STD_NAMES = {"pal": "PAL", "ntsc": "NTSC", "pal_m": "PAL-M", "pal_n": "PAL-N", "secam": "SECAM"}


class AnalogVtxDecoder(EmissionReceiver):
    name = "analog_vtx"
    protocol = "analog FM video (PAL/NTSC)"
    key = "analog_vtx"
    method = "receivers.analog_vtx (research receiver) + FM discriminator"
    owns_band = True                 # FM sidebands / audio subcarrier would otherwise become 'unknown' boxes
    steps = ["~1 us x ~1 MHz time-frequency power; fold every column on the PAL (15625 Hz) and NTSC "
             "(15734.26 Hz) line period over a +-500 ppm grid: the sync tip parks the carrier at one "
             "frequency for 4.7 us per line -> line-locked z-score (works to ~-12 dB in-band SNR)",
             "score >= 10 -> emitter: line rate refined (fold sharpness + drift fit), sync-tip frequency, "
             "occupied band, field timing (lines per field from the vertical sync), standard channel table",
             "FM discriminator on the occupied band: deviation, polarity, audio subcarrier",
             "picture: discriminator output sampled on the measured line clock, active lines stacked "
             "-> <report>_sig<id>_video.png"]

    MIN_RATE_HZ = 20e6               # schema/requirements.json analog_vtx.min_sample_rate_hz: a 6+ MHz FM video

    def run_detector(self, ctx):
        from .detector import detect
        rate = min(ctx.cap.fs, ctx.cap.bandwidth_hz or ctx.cap.fs)
        if rate < self.MIN_RATE_HZ:
            return []                    # cannot hold the signal (SAMPLE_RATE_TOO_LOW warning); a 1 s narrow dwell
                                         # would also cost minutes in the line fold
        return detect(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc, {})

    def in_range(self, lo, hi):
        return _in_ranges(lo, hi, ANALOG_VTX_RANGES)

    def group_key(self, e):
        return e.center_mhz

    def summarise(self, em, sig):
        e = em[0]
        F = e.features
        std = F.get("standard", e.label)
        keys = ("standard", "line_rate_hz", "line_rate_ppm", "line_rate_ppm_sigma", "lines_per_field",
                "field_rate_hz", "scan", "sync_tip_mhz", "center_mhz", "deviation_mhz", "deviation_source",
                "polarity", "audio_subcarrier_mhz", "channel", "channel_candidates", "nominal_mhz",
                "cfo_hz", "bw_mhz", "score", "color", "color_standard", "chroma_fsc_hz", "burst_swing_deg",
                "chroma_coherence", "color_note", "chroma_candidates", "chroma_gate_note", "channel_note",
                "secam_db_hz", "secam_dr_hz", "audio_subcarrier_hz", "audio_dev_rms_khz", "audio_peak_hz",
                "audio_peak_db", "audio_flatness")
        fields = _clean({k: F.get(k) for k in keys if F.get(k) is not None})
        fields["model"] = _model(F, e.center_mhz)          # model level (research, PLAN M8)
        ch = f", channel {F['channel']}" if F.get("channel") else ""
        return result("success", self.method, self.steps, fields,
                      "", label=f"Analog FPV video ({STD_NAMES.get(e.label, STD_NAMES.get(std, 'FM video'))})",
                      confidence=e.confidence, modulation="FM (analog composite video)",
                      evidence=f"line-locked sync tips at {F.get('line_rate_hz', 0):.2f} Hz "
                               f"({F.get('line_rate_ppm', 0):+.0f} ppm), fold score {F.get('score', 0):.1f} "
                               f"(threshold 10){ch}; {F.get('lines_per_field', '?')} lines/field",
                      span=(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6))

    def finish(self, sig, em, ctx):
        from . import detector as M
        try:
            img = M.video_frame(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc, em[0])
        except Exception as e:           # noqa: BLE001 — picture is a bonus, the detection stands
            sig.decoded["fields"]["video_frame_error"] = f"{type(e).__name__}: {e}"[:200]
            return
        if img is not None:
            sig.features["_video_frame"] = img
            sig.decoded["fields"]["video_frame_lines"] = int(img.shape[0])
        sig.decoded["reason"] = ("analog video demodulated to a picture (no CRC exists for analog video; "
                                 "line/field timing locks are the verification)")

    WEAK_DECODERS = {"rc_link", "digital_vtx"}

    def arbitrate(self, signals, ctx):
        """An analog VTX transmits for the whole capture and its FM spectrum (sidebands,
        audio subcarrier, field-blanking transients) reaches beyond its 99.5 % band; weak
        receiver detections (blind OFDM numerology, uncorroborated FSK) inside it are
        pieces of the VTX -- the line period alone looks like a 16 kHz 'OFDM' symbol to the
        blind numerology estimator. Same margin as research drone_survey.explain_by_continuous:
        max(6 MHz, 0.2 x VTX bandwidth). CRC-verified and Herelink (PSS/SSS) signals are kept."""
        vtx = [s for s in signals if s.features.get("decoder") == self.name]
        if not vtx:
            return signals
        keep = []
        for s in signals:
            host = None
            if s.features.get("decoder") in self.WEAK_DECODERS and s.decoded.get("status") != "success" \
                    and "Herelink" not in s.label:
                host = next((v for v in vtx if v.f0 - max(6e6, 0.2 * v.bw) <= s.fc <= v.f1 + max(6e6, 0.2 * v.bw)),
                            None)
            if host is None:
                keep.append(s)
                continue
            host.detection.setdefault("explained_receiver_detections", []).append(
                {"decoder": s.features.get("decoder"), "label": s.label, "f_low_hz": round(s.f0),
                 "f_high_hz": round(s.f1), "note": "inside the analog VTX's FM spectrum: attributed to it"})
        return keep


PLUGINS = [AnalogVtxDecoder()]


def _model(F, center_mhz):
    """Analog VTX model block (research/drone_transmitters/11_model_level_identification.md section 3): the RF shows
    the channel and whether a microphone audio subcarrier is sent; control (SmartAudio / Tramp), power setting and the
    camera standard are not model features, so the model stays unresolved unless a channel plan is model-specific."""
    from spectrum_intel.core.drone_models import model_block
    cands, ev = [], []
    if F.get("channel"):
        ev.append(f"channel {F['channel']}")
    if abs(center_mhz - 5584.0) <= 2.0:
        cands.append("TBS Unify Pro in pit mode (very low power)")
        ev.append("5584 MHz is the TBS Unify Pro pit-mode frequency")
    if any(abs(center_mhz - f) <= 3.0 for f in (1258.0, 1280.0)):
        cands += ["2-channel 1.2 / 1.3 GHz VTX (VTX-1G3SE class)", "9-channel 1.2 / 1.3 GHz VTX"]
        ev.append("1258 / 1280 MHz: the two channels of locked 1.2 / 1.3 GHz VTXs")
    a = F.get("audio_subcarrier_mhz")
    ev.append(f"audio subcarrier {a} MHz: a VTX with a microphone input" if a else "no audio subcarrier")
    return model_block("analog FPV VTX", cands, ev, "inferred",
                       None if cands else "analog VTX models are not separable on air (power, control and camera "
                                          "standard are not visible)")
