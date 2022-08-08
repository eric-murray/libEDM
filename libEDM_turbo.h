#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_convolutional.h>

#include <libEDM_interleaver.h>
#include <libEDM_library.h>

class TurboCodec : public Codec {
public:
    double rate() const {return divide(numUncoded, numCoded);}

    TurboCodec(const size_t  feedbackGenerator,
               const uVector       &parityGenerators,
               const size_t         constraintLength,
               const uMatrix       &interleaverSequences,
               const size_t         numIterations        = 8,
               const Metric         metric               = LOGMAX,
               const bool           adaptiveStop         = false)
        : rscCodec        (constraintLength, parityGenerators, feedbackGenerator, true, true),
          numIterations   (numIterations),
          metric          (metric),
          adaptiveStop    (adaptiveStop),
          Lc              (1.0),
          n               (parityGenerators.size()),
          m               (constraintLength - 1),
          numUncoded      (interleaverSequences[0].size()),
          numCoded        (numUncoded*(1+numCoders*n) + numCoders*m*(1+n)),
          numCoders       (interleaverSequences.size()+1)
          {
              for (size_t i=0; i<interleaverSequences.size(); i++)
              {
                  bInterleavers.push_back(bSequenceInterleaver(interleaverSequences[i]));
                  dInterleavers.push_back(dSequenceInterleaver(interleaverSequences[i]));
              }
          }

    void set_awgn_channel_parameters(const double channelSNR)    {set_scaling_factor(4.0 * channelSNR);}
    void set_iterations             (const size_t numIterations) {this->numIterations = numIterations;}
    void set_scaling_factor         (const double Lc)            {this->Lc = Lc;}

    void    encode (const bVector &input, bVector &output) const {output = encode(input);}
    bVector encode (const bVector &input)                  const;

    void    decode (const dVector &input, bVector &output, const bVector &trueBits = bVector()) {output = decode(input, trueBits);}
    bVector decode (const dVector &input, const bVector &trueBits = bVector());

private:
    const size_t n, m, numCoders;
	const size_t  numUncoded, numCoded;
    const bool   adaptiveStop;
    const Metric metric;

    size_t numIterations;
    double Lc;

    vector<bSequenceInterleaver> bInterleavers;
    vector<dSequenceInterleaver> dInterleavers;

    ConvolutionalCodec rscCodec;

    dVector unsplice (dVector::const_iterator rxSignal, dMatrix &received) const;

};

uVector wcdma_turbo_interleaver_sequence(size_t interleaverSize);
