"""devices -- who is transmitting, what they are doing, and what deserves attention.

Not a receiver: it decodes no IQ. It reads the 802.11 frames other receivers decoded, through the
generic report contract (signal.decoded["fields"]["ppdus"]), and turns them into a device picture and
a list of findings. Keeping it separate means the Wi-Fi receivers stay about demodulation, and this
layer can later take in BLE or anything else that yields addressed frames.

    tracker.tracks(frames)        -> one record per transmitting device
    tracker.link_randomized(devs) -> randomised MACs that are probably one device
    threats.findings(devs, frames, links)
"""
from .tracker import frame, link_randomized, tracks
from .threats import findings

__all__ = ["frame", "tracks", "link_randomized", "findings"]
