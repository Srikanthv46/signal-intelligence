"""Research orchestrators (standalone CLIs that compose several receivers):

    python -m spectrum_intel.research.spectrum_survey cap.sigmf-meta   Wi-Fi + 15 kHz OFDM + DroneID survey
    python -m spectrum_intel.research.drone_survey cap.sigmf-meta      all drone receivers + emitter association

They are not used by the main pipeline (spectrum_intel.cli)."""
