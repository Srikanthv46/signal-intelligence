"""Receivers: one self-contained package per signal family.

Each package keeps its DSP in modules that import only numpy/scipy and
spectrum_intel.core, re-exports its public API from __init__, and has a
plugin.py (the only module that imports the pipeline) exposing PLUGINS for
spectrum_intel.decode.registry. To reuse a receiver elsewhere, import its API
(``from spectrum_intel.receivers.analog_vtx import detect, video_frame``) or copy
the package together with spectrum_intel/core. To drop it from the pipeline,
delete the folder or run with ``--disable <name>``.

Nothing here is imported until asked for.
"""
