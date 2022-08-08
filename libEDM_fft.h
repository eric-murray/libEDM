#pragma once

#include <libEDM_matrix.h>

class FFT {
public:
	static cVector fft (const cVector &input, const bool shift = false) {return fft_select(input, false, shift);}
    static cVector ifft(const cVector &input, const bool shift = false) {return fft_select(input, true,  shift);}

private:
	static cVector fft_select   (const cVector &input, const bool inverse, const bool shift);
	static cVector power2_fft   (const cVector &input, const bool inverse, const bool shift);
	static cVector matlab_fft   (const cVector &input, const bool inverse, const bool shift);
    static size_t  reverse_bits (size_t index, const size_t numBits);
    static cVector fft_shift    (const cVector &input);
};