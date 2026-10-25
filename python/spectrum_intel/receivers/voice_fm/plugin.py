"""Pipeline plugin: narrowband FM voice (walkie-talkies / PMR), demodulated to audio (wraps .detector)."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean, _in_ranges

from .detector import VOICE_RANGES


class VoiceFmDecoder(EmissionReceiver):
    name = "voice_fm"
    protocol = "analog narrowband FM voice"
    key = "voice_fm"
    method = "receivers.voice_fm: FM discriminator + CTCSS / DCS (Golay) + channel plan"
    owns_band = True                 # FM sidebands would otherwise become 'unknown' boxes
    steps = ["1 kHz-bin spectrum of the capture: narrow peaks >= 10 dB over the floor and their "
             "surroundings in 118-174 / 225-520 / 806-870 MHz",
             "per peak, the +-10 kHz channel at 48 kSps: keyed time (5 ms frames >= 10 dB over the floor), "
             "constant envelope, FM discriminator -> carrier, peak deviation (<= 6 kHz), voice-band activity; "
             "2/4-level FSK data and unmodulated carriers are rejected; without a sub-audible code the envelope "
             "must vary no more than noise explains (two-tone beats / clock lines are not FM)",
             "the receiver's own clock lines (pipeline: k x fs, reference-locked) are skipped",
             "sub-audible signalling on the longest keyed segment: DCS 134.4 bit/s NRZ, every 23-bit window "
             "checked against Golay(23,12) and the 100 marker (code + aliases), else a CTCSS tone matched to "
             "the EIA list",
             "clipping products / IQ image of a stronger carrier (k x its LO offset, |k| x its deviation) "
             "dropped and listed on it",
             "channel plan (PMR446, LPD433, FRS) with the receiver ppm when the pipeline has one",
             "audio: channel -> discriminator -> 300-3400 Hz -> 750 us de-emphasis, carrier squelch "
             "-> <name>_voice<MHz>_audio.wav"]

    def run_detector(self, ctx):
        from .detector import detect
        cal = ctx.cache.get("ppm") or ctx.cache.get("sweep_ppm")
        return detect(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc, {"ch": ctx.ch, "ppm": cal["ppm"] if cal else 0.0,
                                                           "exclude_hz": [ln["freq_hz"] for ln in
                                                                          ctx.cache.get("receiver_lines", [])]})

    def in_range(self, lo, hi):
        return _in_ranges(lo, hi, VOICE_RANGES)

    def group_key(self, e):
        return round(e.center_mhz, 3)

    def summarise(self, em, sig):
        e = max(em, key=lambda e: e.confidence)
        F = e.features
        if e.label == "am_voice":
            return self._summarise_am(e, em)
        if e.label in ("aprs", "ax25_packet"):
            return self._summarise_packets(e, em)
        keys = ("center_mhz", "deviation_hz", "channel_bw_khz", "service", "channel", "nominal_mhz", "offset_hz",
                "offset_ppm", "receiver_ppm_applied", "type", "dcs_code", "dcs_aliases", "dcs_standard",
                "dcs_words_verified", "dcs_bit_repeat_agreement", "dcs_deviation_hz", "dcs_note", "ctcss_hz",
                "ctcss_measured_hz", "ctcss_snr_db", "ctcss_standard", "ctcss_deviation_hz", "keyed_s",
                "keyed_fraction", "voice_activity", "voice_rms_hz", "voice_kurtosis", "envelope_cv", "snr_db",
                "segments_s", "suppressed_artifacts")
        fields = _clean({k: F.get(k) for k in keys if F.get(k) is not None})
        fields["audio"] = "nbfm"                       # the viewer offers 'Listen' on this
        code = F.get("dcs_code") or (f"CTCSS {F['ctcss_hz']} Hz" if F.get("ctcss_hz") else None)
        where = f"{F['service']} ch {F['channel']}" if F.get("service") else f"{e.center_mhz:.4f} MHz"
        verified = bool(F.get("dcs_code") and F.get("dcs_words_verified", 0) >= 2) or bool(F.get("ctcss_hz"))
        speech = F.get("voice_activity", 0) >= 0.1
        label = f"Analog FM voice ({where}{', ' + code if code else ''})" if (code or speech) \
            else f"Analog FM voice-like ({where})"
        reason = ("FM voice demodulated to audio; " +
                  ("sub-audible squelch code verified (DCS Golay / standard CTCSS tone)" if verified else
                   "no sub-audible code to verify") +
                  ("" if speech else "; little speech-band activity in this capture (carrier keyed, quiet mic?)"))
        return result("success" if verified else "partial", self.method, self.steps, fields, reason,
                      label=label, confidence=e.confidence, modulation="FM (analog narrowband voice)",
                      evidence="; ".join(F.get("evidence", [])),
                      span=(min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6))

    def _summarise_am(self, e, em):
        """AM voice (research, signal gaps S13)."""
        F = e.features
        keys = ("center_mhz", "modulation_rms", "channel_bw_khz", "service", "channel_mhz", "offset_hz", "emergency",
                "keyed_s", "keyed_fraction", "voice_activity", "phase_rms_hz", "envelope_cv", "snr_db", "segments_s",
                "suppressed_artifacts")
        fields = _clean({k: F.get(k) for k in keys if F.get(k) is not None})
        fields["audio"] = "am"                         # the viewer offers 'Listen' on this
        where = F.get("emergency") or (f"{F['service']}, {F['channel_mhz']} MHz" if F.get("service")
                                       else f"{e.center_mhz:.4f} MHz")
        label = f"Analog AM voice ({where})"
        return result("partial", "receivers.voice_fm: AM envelope detector (research)",
                      [self.steps[0], "per peak, the +-10 kHz channel at 48 kSps: keyed time; AM: the envelope's "
                       "300-3400 Hz rms relative to the carrier (>= 4 % and 3x what noise gives) in >= 10 % of 20 ms "
                       "frames, while the discriminator's speech-band rms stays at the noise level (phase quiet)",
                       "aeronautical band and 25 kHz channel (VHF civil 118-137 MHz, UHF military 225-400 MHz); "
                       "121.5 / 243.0 MHz emergency",
                       "audio: channel -> envelope -> 300-3400 Hz, gain from the carrier level, carrier squelch "
                       "-> <name>_voice<MHz>_audio.wav"], fields,
                      "AM voice demodulated to audio; AM carries no code to verify",
                      label=label, confidence=e.confidence, modulation="AM (analog voice, full carrier)",
                      evidence="; ".join(F.get("evidence", [])),
                      span=(min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6))

    def _summarise_packets(self, e, em):
        """AX.25 / APRS packet radio (research, signal gaps S19)."""
        F = e.features
        keys = ("center_mhz", "deviation_hz", "channel_bw_khz", "keyed_s", "snr_db", "stations", "packets",
                "segments_s")
        fields = _clean({k: F.get(k) for k in keys if F.get(k) is not None})
        what = "APRS" if e.label == "aprs" else "AX.25 packet radio"
        st = F["stations"]
        return result("success", "receivers.voice_fm: Bell 202 AFSK + AX.25 HDLC (research)",
                      [self.steps[0], "per peak, the +-10 kHz channel at 48 kSps: FM discriminator",
                       "Bell 202 AFSK: mark 1200 / space 2200 Hz tone energy over one bit, 8 bit phases, NRZI "
                       "decoding, HDLC frames with a valid CRC-16/X.25, AX.25 UI addresses, APRS position / status / "
                       "message fields"], fields, "every packet verified by its CRC-16/X.25",
                      label=f"{what} ({', '.join(st[:3])}{', ...' if len(st) > 3 else ''})", confidence=e.confidence,
                      modulation="AFSK 1200 bit/s (Bell 202) on narrowband FM", evidence="; ".join(F["evidence"]),
                      span=(min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6))

    def finish(self, sig, em, ctx):
        if sig.decoded["fields"].get("packets") is not None:
            return                                   # packet radio: decoded data, no audio
        from .audio import demodulate
        try:
            a = demodulate(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc, sig.decoded["fields"])
        except Exception as e:           # noqa: BLE001 — audio is a bonus, the detection stands
            sig.decoded["fields"]["audio_error"] = f"{type(e).__name__}: {e}"[:200]
            return
        if a is not None and len(a.samples):
            sig.features["_audio"] = (a.samples, a.fs)
            sig.decoded["fields"]["audio_s"] = round(len(a.samples) / a.fs, 3)
            sig.decoded["fields"]["audio_squelch_open"] = round(a.meta["squelch_open"], 3)


PLUGINS = [VoiceFmDecoder()]
