#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_fft.h>
#include <libEDM_modulators.h>

// class Modulator

cVector Modulator::modulate(const bVector &bits) const
{
    cVector output;
    modulate(bits, output);
    return output;
}

dVector Modulator::demodulate(const cVector &received) const
{
    dVector output;
    demodulate(received, output);
    return output;
}

// class BPSK

void BPSK::modulate(const bVector &bits, dVector &output) const
{
    output.resize(bits.size());
    for (size_t i=0; i<bits.size(); i++)
        output[i] = (bits[i] == false ? 1.0 : -1.0);
}

void BPSK::modulate(const bVector &bits, cVector &output) const
{
    output.resize(bits.size());
    for (size_t i=0; i<bits.size(); i++)
        output[i] = (bits[i] == false ? complex<double>(1.0, 0.0) : complex<double>(-1.0, 0.0));
}

void BPSK::demodulate(const dVector &received, bVector &output) const
{
    output.clear();
    for (size_t i=0; i<received.size(); i++)
        output.push_back((received[i] > 0.0) ? false : true);
}

void BPSK::demodulate(const cVector &received, bVector &output) const
{
    output.clear();
    for (size_t i=0; i<received.size(); i++)
        output.push_back((received[i].real() > 0.0) ? false : true);
}

void BPSK::demodulate(const cVector &received, dVector &output) const
{
    output.clear();
    for (size_t i=0; i<received.size(); i++)
        output.push_back(received[i].real());
}

// class QPSK

void QPSK::modulate(const bVector &bits, cVector &output) const
{
    // check bits contains an even number of bits
    assert( bits.size() % 2 == 0 );

    output.resize(bits.size()/2);
    for (size_t i=0; i<output.size(); i++)
    {
        double real_part = (bits[2*i]   == false) ? M_SQRT1_2 : -M_SQRT1_2;
        double imag_part = (bits[2*i+1] == false) ? M_SQRT1_2 : -M_SQRT1_2;
        output[i] = complex<double>(real_part, imag_part);
    }
}

void QPSK::demodulate(const cVector &received, bVector &output) const
{
    output.clear();
    for (size_t i=0; i<received.size(); i++)
    {
        output.push_back((received[i].real() > 0.0) ? false : true);
        output.push_back((received[i].imag() > 0.0) ? false : true);
    }
}

void QPSK::demodulate(const cVector &received, dVector &output) const
{
    output.clear();
    for (size_t i=0; i<received.size(); i++)
    {
        output.push_back(received[i].real());
        output.push_back(received[i].imag());
    }
}

// class OFDM
void OFDM::modulate(const bVector &bits, cVector &output) const
{
    // check bits contains 2*numFFT bits
    assert( bits.size() == 2 * numFFT );

    // map data bits to complex vector
    cVector input;
    for (size_t i=0; i<numFFT; i++)
        input.push_back(complex<double>(bits[2*i], bits[2*i+1]));

    output = FFT::ifft(input);

    // add cyclic prefix
    output.ins(0, output.right(numCyclicPrefix));
}

void OFDM::demodulate(const cVector &received, bVector &output) const
{
    dVector realOutput;
    demodulate(received, realOutput);

    for (size_t i=0; i<2*numFFT; i++)
        output.push_back((realOutput[i] > 0.0) ? false : true);
}

void OFDM::demodulate(const cVector &received, dVector &output) const
{
    output.clear();

    cVector complexOutput = FFT::fft(received.right(numFFT));

    for (size_t i=0; i<numFFT; i++)
    {
        output.push_back(complexOutput[i].real());
        output.push_back(complexOutput[i].imag());
    }
}




