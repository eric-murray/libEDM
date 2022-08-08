#pragma once

#include <libEDM_matrix.h>

class Modulator {
public:
    Modulator (const size_t bits_per_symbol) : _M(bits_per_symbol) {}

    dVector demodulate (const cVector &received) const;
    cVector modulate   (const bVector &bits)     const;

    size_t bits_per_symbol () const { return _M; }

    virtual void demodulate (const cVector &received, dVector &output) const = 0;
    virtual void modulate   (const bVector &bits, cVector &output)     const = 0;

private:
    const size_t _M;
};

class BPSK : public Modulator {
public:
    BPSK () : Modulator(1) {}

    void demodulate (const dVector &received, bVector &output) const;
    void demodulate (const dVector &received, dVector &output) const { output = received; }
    void demodulate (const cVector &received, bVector &output) const;
    void demodulate (const cVector &received, dVector &output) const;

    void modulate   (const bVector &bits, dVector &output) const;
    void modulate   (const bVector &bits, cVector &output) const;

    dVector demodulate (const cVector &received) const { return Modulator::demodulate(received); }
    cVector modulate   (const bVector &bits)     const { return Modulator::modulate(bits); }
};

class QPSK : public Modulator {
public:
    QPSK () : Modulator(2) {}

    void demodulate (const cVector &received, bVector &output) const;
    void demodulate (const cVector &received, dVector &output) const;

    void modulate   (const bVector &bits, cVector &output) const;

    dVector demodulate (const cVector &received) const { return Modulator::demodulate(received); }
    cVector modulate   (const bVector &bits)     const { return Modulator::modulate(bits); }
};

class OFDM : public Modulator {
public:
    const size_t numFFT, numCyclicPrefix;

    OFDM (size_t numFFT, size_t numCyclicPrefix) : Modulator(2*numFFT), numFFT(numFFT), numCyclicPrefix(numCyclicPrefix) {}
      
    void demodulate (const cVector &received, bVector &output) const;
    void demodulate (const cVector &received, dVector &output) const;

    void modulate   (const bVector &bits, cVector &output) const;

    dVector demodulate (const cVector &received) const { return Modulator::demodulate(received); }
    cVector modulate   (const bVector &bits)     const { return Modulator::modulate(bits); }
};