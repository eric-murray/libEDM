#pragma once

#include <libEDM_types.h>
#include <libEDM_filter.h>
#include <libEDM_library.h>
#include <libEDM_matrix.h>
#include <libEDM_random.h>

class Channel {
public:
    const double awgnStdDev;

protected:
    Channel (const double awgnVariance = 0.0) : awgnStdDev(sqrt(awgnVariance)) {}
};

class AWGNChannel : public Channel {
public:
    AWGNChannel (const double awgnVariance = 0.0) : Channel(awgnVariance) {}

    cVector operator()(const cVector &input);
    dVector operator()(const dVector &input);
};

class RayleighFadingProcess {
public:
    cVector operator()(const cVector &input, cVector &channel);

    cVector fade      (const cVector &input, cVector &channel) {return operator()(input, channel);}
    cVector operator()(const cVector &input)                   {cVector channel; return operator()(input, channel);}
    cVector fade      (const cVector &input)                   {cVector channel; return operator()(input, channel);}

private:
    virtual cVector compute_channel (const size_t numSamples) = 0;
};

class SinusoidRFP : public RayleighFadingProcess {
public:
    SinusoidRFP (const double normalisedDoppler, const double weight, const size_t numSinusoids = 256);

private:
    size_t   _numSinusoids;
    double  _amplitude;
    dVector _frequencies, _phaseOffsets;
    size_t   _index;

    cVector compute_channel (const size_t numSamples);
};

class FilteredRFP : public RayleighFadingProcess {
public:
    FilteredRFP (const double normalisedDoppler, const double weight, const size_t numTaps = 1023) : _fadingFilter(coefficients(numTaps, normalisedDoppler, weight)) {}

private:
    Filter<complex<double>,double> _fadingFilter;

    dVector coefficients    (const size_t numTaps, const double normalisedDoppler, const double weight) const;
    cVector compute_channel (const size_t numSamples);
};

template <class RFP_Type = FilteredRFP>
class MultipathRayleighChannel  : public AWGNChannel {
public:
    enum Model {Rayleigh, Case1, Case3};

    MultipathRayleighChannel (const double awgnVariance, const double doppler, const double samplePeriod, const uVector tapDelays, const dVector tapWeights, const bool normaliseTapWeights = true);
    MultipathRayleighChannel (const double awgnVariance, const double sampleInterval, const double speed_kmh, const double wavelength, Model model, const bool normaliseTapWeights = true);

    cVector operator()(const cVector &input);
    cVector operator()(const cVector &input, cMatrix &channel);

private:
    const double _sampleInterval;       // sample period in seconds

    cVector _residualOutput;

    double  _normalisedDoppler;
    uVector _tapDelays;
    dVector _tapWeights, _normalisedLinearTapWeights;

    vector<RFP_Type*> rayleighFaders;
};

class OFDMChannel {
public:
    enum Model {WINNER_SISO_URBAN_MACRO, WINNER_SISO_URBAN_MICRO};

    const double doppler;

    double time() const {return _time;}

    double  frequency_response (const double frequency_Hz) const;
    dVector frequency_response (const double startFrequency_Hz, const double stepFrequency_Hz, const size_t numFrequencies) const;
    dVector frequency_response (const double frequencyInterval_Hz, const double timeInterval_s) const;

    OFDMChannel (const Model model, const double speed_kmh, const double frequency_Hz);

    void set_time     (const double time);
    void advance_time (const double deltaTime) {set_time(_time + deltaTime);}

private:
    class RayleighFadingProcess {
    public:
        const double delay;

        complex<double> gain () {return _gain;}

        RayleighFadingProcess (const double doppler, const double weight, const double delay, const size_t numSinusoids = 20);

        void update (const double time);

    private:
        size_t          _numSinusoids;
        double          _amplitude;
        dVector         _frequencies, _phaseOffsets;
        double          _time;
        complex<double> _gain;
    };

    double  _time;
    dVector _tapDelays;
    dVector _tapWeights;

    vector<RayleighFadingProcess*> _taps;

};

//
// Templated class implementations
//

//
// class MultipathRayleighChannel
//

template <class RFP_Type>
MultipathRayleighChannel<RFP_Type>::MultipathRayleighChannel (double awgnVariance, double doppler, double sampleInterval, uVector tapDelays, dVector tapWeights, bool normaliseTapWeights)
: AWGNChannel(awgnVariance), _normalisedDoppler(doppler * sampleInterval), _sampleInterval(sampleInterval), _tapDelays(tapDelays), _tapWeights(tapWeights)
{
    for (size_t i=0; i<_tapWeights.size(); i++)
        _normalisedLinearTapWeights.push_back(pow(10.0, 0.05*_tapWeights[i]));

    if ( normaliseTapWeights )
        _normalisedLinearTapWeights / _normalisedLinearTapWeights.sum_of_squares();

    for (size_t i=0; i<_tapDelays.size(); i++)
        rayleighFaders.push_back(RFP_Type(_normalisedDoppler, _normalisedLinearTapWeights[i]));
}

template <class RFP_Type>
MultipathRayleighChannel<RFP_Type>::MultipathRayleighChannel (double awgnVariance, double sampleInterval, double speed_kmh, double wavelength, Model model, bool normaliseTapWeights)
: AWGNChannel(awgnVariance), _sampleInterval(sampleInterval)
{
    _normalisedDoppler = speed_kmh / (3.6 * wavelength) * sampleInterval;

    switch ( model )
    {
    case Rayleigh:
        _tapDelays .resize(1);
        _tapWeights.resize(1);

        // tap 0
        _tapDelays [0] = 0;
        _tapWeights[0] = 0.0;

        break;

    case Case1:
        _tapDelays .resize(2);
        _tapWeights.resize(2);

        // tap 0
        _tapDelays [0] = 0;
        _tapWeights[0] = 0.0;

        // tap 1
        _tapDelays [1] = round(976E-9 / sampleInterval);
        _tapWeights[1] = -10.0;

        break;

    case Case3:
        _tapDelays .resize(4);
        _tapWeights.resize(4);

        // tap 0
        _tapDelays [0] = 0;
        _tapWeights[0] = 0.0;

        // tap 1
        _tapDelays [1] = round(260E-9 / sampleInterval);
        _tapWeights[1] = -3.0;

        // tap 2
        _tapDelays [2] = round(521E-9 / sampleInterval);
        _tapWeights[2] = -6.0;

        // tap 3
        _tapDelays [3] = round(781E-9 / sampleInterval);
        _tapWeights[3] = -9.0;

        break;
    }

    for (size_t i=0; i<_tapWeights.size(); i++)
        _normalisedLinearTapWeights.push_back(pow(10.0, 0.05*_tapWeights[i]));

    if ( normaliseTapWeights )
        _normalisedLinearTapWeights / _normalisedLinearTapWeights.sum_of_squares();

    for (size_t i=0; i<_tapDelays.size(); i++)
        rayleighFaders.push_back(new RFP_Type(_normalisedDoppler, _normalisedLinearTapWeights[i]));
}

template <class RFP_Type>
cVector MultipathRayleighChannel<RFP_Type>::operator()(const cVector &input, cMatrix &channel)
{
    const size_t length   = input.size();
    const size_t maxDelay = _tapDelays.back();

    // copy residual output from previous filering to current output
    cVector output = _residualOutput;
    output.resize(length);

    // reset _residualOutput
    _residualOutput.assign(maxDelay, complex<double>(0.0, 0.0));

    channel.set_size(_tapDelays.size(), 0, 0.0);
    for (size_t i=0; i<_tapDelays.size(); i++)
    {
        // compute contribution of tap i to output
        cVector tapContribution(length + maxDelay);
        tapContribution.replace_mid(_tapDelays[i], rayleighFaders[i]->fade(input, channel[i]));

        // update output and residualOutput
        output          += tapContribution.left (length);
        _residualOutput += tapContribution.right(maxDelay);
    }

    return AWGNChannel::operator ()(output);
}

template <class RFP_Type>
cVector MultipathRayleighChannel<RFP_Type>::operator()(const cVector &input)
{
    cMatrix channel;
    return operator()(input, channel);
}