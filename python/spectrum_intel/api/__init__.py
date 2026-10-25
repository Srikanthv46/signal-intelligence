"""Production API of spectrum_intel: IQ blocks in, versioned JSON responses out.

    from spectrum_intel.api import Engine, IQBlock, StreamAnalyzer, analyze_block

    resp = analyze_block(iq, sample_rate_hz=61.44e6, center_frequency_hz=5.8e9, bandwidth_hz=56e6)

The response schema is schema/response.schema.json (shared with the C++ library);
requirements() gives the sample rate and block length each receiver needs.
"""
from .contract import SCHEMA_VERSION, requirements, technology
from .engine import Engine, IQBlock, analyze_block, to_complex64, validate
from .stream import StreamAnalyzer

__all__ = ["SCHEMA_VERSION", "Engine", "IQBlock", "StreamAnalyzer", "analyze_block", "requirements", "technology",
           "to_complex64", "validate"]
