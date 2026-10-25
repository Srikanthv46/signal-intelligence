"""Analog FPV video (FM PAL/NTSC composite) receiver.

    detect(x, fs, fc, ctx=None) -> list[Emission]      family "analog_vtx", label pal | ntsc | analog_fm_video
    video_frame(x, fs, fc, emission) -> uint8 [lines, 720] picture (or None)
    FieldDemodulator.from_fields(fields).process(iq, fs, fc) -> VideoFrame   streaming video, one field per block
    synth                                               synthetic VTX generator (tests)

Depends on spectrum_intel.core only. See README.md.
"""
from .detector import CHANNELS, THRESHOLD, detect, match_channel, video_frame  # noqa: F401
from .video import FieldDemodulator, VideoFrame  # noqa: F401
