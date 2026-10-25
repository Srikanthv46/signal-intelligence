"""Wi-Fi receivers.

    ofdm_rx.WifiDetector(fs, fc).run(x) -> packets     802.11a/g/n/ac/ax legacy preamble, L-SIG, HT/VHT-SIG
    payload.decode_nonht_payload(...)                  non-HT MPDU, FCS, MAC header, SSID
    payload.decode_ht_payload(...)                     802.11n single stream: BCC/LDPC, with or without STBC
    ldpc.params / decode_stream / encode_stream        802.11 LDPC codes and PPDU codeword layout
    ofdm_scaled.scan(x, fs, fc) -> [Link]              802.11 at another clock / frequency / mirrored, any band
    links.attribute(records) -> transmitters           which station sent each PPDU (MAC, responses, seed counters)
    dsss.scan(IQContext(iq, fs, fc)) -> list           802.11b DSSS/CCK: PLCP CRC-16, MPDU FCS, MAC dissection
    mac80211.parse(mpdu) -> dict                       802.11 MAC dissector (IEs, security, WPS, EAPOL, ...)
    dsss_tx.ppdu_chips / beacon / scene                reference 802.11b transmitter + channel model (tests)
    ofdm_tx.ppdu / ampdu / timescale / scene           reference 802.11a/g/n transmitter (LDPC, STBC, any clock)

Depends on spectrum_intel.core only. See README.md.
"""
from .ofdm_rx import WifiDetector  # noqa: F401
from .payload import decode_nonht_payload  # noqa: F401
