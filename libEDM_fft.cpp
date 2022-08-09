#define _USE_MATH_DEFINES

#include <complex>
#include <fstream>
#include <iostream>

#include <libEDM_fft.h>
#include <libEDM_library.h>
#include <libEDM_matlab.h>
#include <libEDM_matrix.h>

using std::cerr;
using std::polar;

cVector FFT::fft_select (const cVector &input, const bool inverse, const bool shift)
{
    if ( pow(2.0, ceil(log2(input.size()))) == input.size() )
        // FFT length is a power of 2
        return power2_fft(input, inverse, shift);
    else
        return matlab_fft(input, inverse, shift);
}

cVector FFT::matlab_fft (const cVector &input, const bool inverse, const bool shift)
{
    // check that MATLAB engine has been started
    if ( !MATLAB::engine )
    {
        MATLAB::openEngine();
    }

    // create mxArray and copy input values to it
    mxArray *inputArray = convert_mxType(input);

    // put inputArray into MATLAB workspace
    engPutVariable(MATLAB::engine, "input", inputArray);

    // compute fft or ifft of inputArray, and store in variable output
    if ( inverse )
        if ( shift )
   	        engEvalString(MATLAB::engine, "output = ifft(fftshift(input))");
        else
   	        engEvalString(MATLAB::engine, "output = ifft(input)");
    else
        if ( shift )
   	        engEvalString(MATLAB::engine, "output = fftshift(fft(input))");
        else
   	        engEvalString(MATLAB::engine, "output = fft(input)");

    // create mxArray and copy fftarray data to it
    mxArray *outputArray = engGetVariable(MATLAB::engine, "output");

    // create output vector and copy mxArray data to it
    cVector output = convert_mxType(outputArray);

    // destroy mxArrays to release memory
    mxDestroyArray(inputArray);
    mxDestroyArray(outputArray);

    return output;
}

cVector FFT::power2_fft (const cVector &input, const bool inverse, const bool shift)
{
    // Pad input data out to higher power of two
    const size_t numBits    = static_cast<size_t>(ceil(log2(input.size())));
    const size_t numSamples = static_cast<size_t>(pow(2.0, numBits));

    // Do simultaneous data copy and bit-reversal ordering into outputs
    // If inverse transform, and shift requested, do that as well
    cVector output(numSamples);
    for (size_t i=0; i < input.size(); i++)
    {
        size_t j = reverse_bits(i, numBits);
        if ( !inverse || !shift )
		    output[j] = input[i];
        else
            if ( i < input.size()/2 )
                output[j] = input[i + input.size()/2];
            else
                output[j] = input[i - input.size()/2];
    }

    size_t blockEnd = 1;
    for (size_t blockSize = 2; blockSize < 2*numSamples; blockSize *= 2)
    {
		const double theta = (inverse ? 2.0 : -2.0) * M_PI / blockSize;
        complex<double> wp = polar(1.0, theta);

        complex<double> w(1.0, 0.0);
		for (size_t m=0; m<blockEnd; m++)
        {
			for (size_t i=m; i<numSamples; i+=blockSize)
            {
				const size_t j = i + blockEnd;
                const complex<double> temp = w * output[j];

                output[j]  = output[i] - temp;
                output[i] += temp;
			}
            w *= wp;
		}
        blockEnd = blockSize;
	}

    // Need to normalize if inverse transform
    if ( inverse )
		output /= numSamples;
    else
        if ( shift )
            output = fft_shift(output);

    return output;
}

cVector FFT::fft_shift (const cVector &input)
{
    size_t halfLength = input.size()/2;

    cVector output;
    output.reserve(input.size());
    output.ins(0, input.right(halfLength));
    output.ins(halfLength, input.left(halfLength));

    return output;
}

size_t FFT::reverse_bits (size_t index, const size_t numBits)
{
    size_t reverse = 0;
    for (size_t i=0; i < numBits; i++ )
    {
        reverse = (reverse << 1) | (index & 1);
        index >>= 1;
    }

    return reverse;
}
