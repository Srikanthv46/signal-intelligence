// Short-time power spectra on absolute frequency / time axes (python: detect/stft.py).
#pragma once

#include "core/spectrogram.hpp"

namespace si::detect {

// Hann-windowed, non-padded frames every hop samples (hop = nfft when 0); P = |fftshift(FFT)|^2 / sum(w^2),
// float32, so white noise of variance s2 gives E[P] = s2. f: absolute bin centres; t: frame starts.
Spectrogram stft(const cfvec& x, double fs, double fc, size_t nfft, size_t hop = 0);

}  // namespace si::detect
