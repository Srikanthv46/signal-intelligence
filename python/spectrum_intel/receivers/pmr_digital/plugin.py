"""Pipeline plugin: digital land-mobile radio (TETRA, DMR, NXDN), wraps .detector."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean, _in_ranges

from .detector import FSK_RANGES, TETRA_RANGES

_TETRA_KEYS = ("carrier_hz", "mcc", "mnc", "colour_code", "system_code", "sharing_mode", "ts_reserved_frames",
               "u_plane_dtx", "frame18_extension", "cell_service_level", "late_entry", "neighbour_cell_broadcast",
               "n_bsch_crc_ok", "n_sync_bursts", "n_normal_training", "slot_consistent_hits", "bsch",
               "bsch_identities", "symbol_rate_line", "phase_error_rad", "snr_db")
_DMR_KEYS = ("source", "colour_code", "colour_code_votes", "n_bursts", "sync_types", "voice_bursts", "data_types",
             "crc_checked", "crc_ok", "csbk", "symbol_rate", "symbol_rate_line", "deviation_hz_p99", "snr_db", "bursts")
_NXDN_KEYS = ("symbol_rate", "channel_bw_khz", "n_frame_sync", "frame_consistent_hits", "symbol_rate_line",
              "deviation_hz_p99", "snr_db")


class PmrDigitalDecoder(EmissionReceiver):
    name = "pmr_digital"
    protocol = "digital land-mobile radio (TETRA, DMR, NXDN)"
    key = "pmr_digital"
    method = "receivers.pmr_digital: TETRA pi/4-DQPSK + BSCH; 4FSK discriminator + DMR sync/Golay/BPTC, NXDN sync"
    owns_band = True
    steps = ["candidates: narrow carriers (<= 30 kHz, >= 8 dB) in 380-400/410-430/450-470/870-876/915-921 MHz "
             "(TETRA) and 136-174/300-520/806-870 MHz (4FSK), 1 kHz bins",
             "TETRA: 25 kHz channel at 144 kSps, 18 ksym/s envelope line; pi/4-DQPSK (RRC 0.35) with the "
             "4th-power carrier estimate and its quarter-turn ambiguity resolved by the training sequences; "
             "synchronisation (38-bit) and normal (22-bit) training sequences on the 510-bit slot grid; "
             "BSCH: descramble (e=0), de-interleave (120, 11), RCPC 2/3 -> rate-1/4 K=5 soft Viterbi, CRC-CCITT "
             "-> SYNC PDU (MCC, MNC, colour code, TN/FN/MN)",
             "4FSK: 12.5 kHz channel at 48 kSps, FM discriminator, 2400/4800 sym/s line; DMR: SYNC (9 patterns) "
             "-> per-burst RRC timing/levels from the sync -> slot type Golay(20,8) (colour code, data type) "
             "-> BPTC(196,96) Hamming(13,9)/(15,11) -> CRC-CCITT with the data-type mask (CSBK opcode when valid); "
             "NXDN: 20-bit frame sync on the 192-symbol frame"]

    def run_detector(self, ctx):
        from .detector import detect
        det = getattr(ctx, "det", None)
        return detect(ctx.ch, {"band": getattr(det, "usable_hz", None),
                               "exclude_hz": [ln["freq_hz"] for ln in ctx.cache.get("receiver_lines", [])]})

    def in_range(self, lo, hi):
        return _in_ranges(lo, hi, TETRA_RANGES + FSK_RANGES)

    def group_key(self, e):
        return (e.family, round(e.center_mhz, 3))

    def keep(self, e):
        return True

    def summarise(self, em, sig):
        e = max(em, key=lambda e: e.confidence)
        F = e.features
        span = (min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6)
        where = f"{e.center_mhz:.4f} MHz"
        if e.family == "tetra":
            fields = _clean({k: F.get(k) for k in _TETRA_KEYS if F.get(k) is not None})
            if e.label == "tetra_downlink_bsch":
                lab = f"TETRA downlink MCC {F['mcc']} MNC {F['mnc']} colour code {F['colour_code']} ({where})"
                return result("success", self.method, self.steps, fields,
                              "BSCH (SYNC PDU) CRC-verified; MAC/MLE broadcast (SYSINFO) not decoded", label=lab,
                              confidence=e.confidence, modulation="pi/4-DQPSK 18 ksym/s (TETRA)",
                              evidence=f"{F['n_bsch_crc_ok']} of {F['n_sync_bursts']} synchronisation bursts' BSCH "
                                       f"pass the CRC; training sequences on one slot phase "
                                       f"({F['slot_consistent_hits']}); 18 ksym/s line x{F['symbol_rate_line']}",
                              span=span)
            if e.label == "tetra_downlink":
                return result("partial", self.method, self.steps, fields,
                              "TETRA burst structure (training sequences) found; no synchronisation burst decoded",
                              label=f"TETRA downlink ({where})", confidence=e.confidence,
                              modulation="pi/4-DQPSK 18 ksym/s (TETRA)",
                              evidence=f"{F['slot_consistent_hits']} training sequences on one 510-bit slot phase; "
                                       f"18 ksym/s line x{F['symbol_rate_line']}", span=span)
            return result("partial", self.method, self.steps, fields,
                          "18 ksym/s pi/4-DQPSK-like carrier, no training sequence in this block",
                          label=f"TETRA-like carrier ({where})", confidence=e.confidence,
                          modulation="pi/4-DQPSK-like 18 ksym/s",
                          evidence=f"18 ksym/s line of the matched-filtered envelope x{F['symbol_rate_line']}",
                          span=span)
        if e.family == "dmr":
            fields = _clean({k: F.get(k) for k in _DMR_KEYS if F.get(k) is not None})
            cc = f", colour code {F['colour_code']}" if F.get("colour_code") is not None else ""
            dts = F.get("data_types") or {}
            what = ("voice call" if F.get("voice_bursts") else
                    ("control channel (CSBK)" if dts and max(dts, key=dts.get) == "CSBK" else ""))
            lab = f"DMR {F['source']}{cc}{', ' + what if what else ''} ({where})"
            ok = bool(F.get("csbk"))
            reason = ("CSBK CRC-verified" if ok else
                      "SYNC, slot type (Golay) and BPTC (Hamming) verified; " +
                      (f"CSBK/header CRC failed {F['crc_checked'] - F.get('crc_ok', 0)} of {F['crc_checked']} "
                       "(non-standard or vendor control channel?)" if F.get("crc_checked") else
                       "no CRC-protected block in this capture"))
            return result("success" if ok else "partial", self.method, self.steps, fields, reason, label=lab,
                          confidence=e.confidence, modulation="4FSK 4800 sym/s (DMR)",
                          evidence=f"{F['n_bursts']} DMR SYNC word(s) {F['sync_types']}" +
                                   (f"; colour code {F['colour_code']} in {F['colour_code_votes']} slot types "
                                    f"(Golay)" if F.get("colour_code") is not None else ""), span=span)
        fields = _clean({k: F.get(k) for k in _NXDN_KEYS if F.get(k) is not None})
        return result("partial", self.method, self.steps, fields,
                      "NXDN frame sync verified; LICH/SACCH not decoded",
                      label=f"NXDN {F['symbol_rate']:g} sym/s ({where})", confidence=e.confidence,
                      modulation=f"4FSK {F['symbol_rate']:g} sym/s (NXDN)",
                      evidence=f"{F['n_frame_sync']} NXDN frame sync word(s), {F['frame_consistent_hits']} on the "
                               f"192-symbol frame grid; symbol-rate line x{F['symbol_rate_line']}", span=span)


    def arbitrate(self, signals, ctx):
        """A carrier verified here (sync words / CRC) is digital: drop another receiver's analog reading of it
        (voice_fm: '4FSK' passes as FM voice in short blocks)."""
        mine = [s for s in signals if s.features.get("decoder") == self.name]
        if not mine:
            return signals
        return [s for s in signals if s.features.get("decoder") != "voice_fm" or
                not any(abs(s.fc - m.fc) < 4e3 for m in mine)]


PLUGINS = [PmrDigitalDecoder()]
