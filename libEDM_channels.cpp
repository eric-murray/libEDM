#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <fstream>

#include <itpp/base/bessel/bessel_internal.h>

#include <libEDM_channels.h>
#include <libEDM_fft.h>
#include <libEDM_matrix.h>

using std::abs;
using std::cos;
using std::norm;
using std::polar;

//
// class AWGNChannel
//

dVector AWGNChannel::operator()(const dVector &input)
{
    dVector output;
    for (size_t i=0; i<input.size(); i++)
        output.push_back(input[i] + globalRandom.gaussian(0.0, awgnStdDev));
    return output;
}

cVector AWGNChannel::operator()(const cVector &input)
{
    cVector output;
    for (size_t i=0; i<input.size(); i++)
        output.push_back(input[i] + globalRandom.complex_gaussian(0.0, awgnStdDev));
    return output;
}

//
// class RayleighFadingProcess
//

cVector RayleighFadingProcess::operator() (const cVector &input, cVector &channel)
{
    channel = compute_channel(input.size());

    cVector output = input * channel;
    return output;
}

//
// class SinusoidRFP
//

SinusoidRFP::SinusoidRFP (const double normalisedDoppler, const double weight, const size_t numSinusoids)
: _numSinusoids(numSinusoids), _index(0), _amplitude(weight / sqrt(static_cast<double>(_numSinusoids)))
{
    _frequencies.resize(numSinusoids, 0.0);
    for (size_t i=0; i<_numSinusoids; i++)
        _frequencies[i] = -2.0 * M_PI * normalisedDoppler * cos(globalRandom.angle());

    _phaseOffsets.resize(numSinusoids, 0.0);
    for (size_t i=0; i<_numSinusoids; i++)
        _phaseOffsets[i] = globalRandom.angle();
}

cVector SinusoidRFP::compute_channel (const size_t numSamples)
{
    cVector channel(numSamples, 0.0);
    for (size_t i=0; i<_numSinusoids; i++)
        for (size_t n=0; n<numSamples; n++)
        {
            const double phase = _frequencies[i]*(_index+n) - _phaseOffsets[i];
            channel[n] += polar(_amplitude, phase);
        }

    _index += numSamples;

    return channel;
}

//
// class FilteredRFP
//

cVector FilteredRFP::compute_channel (const size_t numSamples)
{
    // compute complex gaussian noise samples
    cVector gaussianNoise;
    for (size_t i=0; i<numSamples; i++)
        gaussianNoise.push_back(globalRandom.complex_gaussian());

    // compute channel by filtering gaussian noise and return
    return _fadingFilter(gaussianNoise);
}

dVector FilteredRFP::coefficients (const size_t numTaps, const double normalisedDoppler, const double weight) const
{
    dVector coefficients(numTaps);

    const size_t centreTap = (numTaps - 1) / 2;
    const double constant  = pow(normalisedDoppler * M_1_PI, 0.25) * gam(0.75);

    for (size_t n=0; n<numTaps; n++)
    {
        if ( n != centreTap )
        {
            const double temp = fabs(static_cast<double>(n) - centreTap);
            coefficients[n] = constant * pow(temp, -0.25) * jv(0.25, 2.0 * M_PI * normalisedDoppler * temp);
        }
        else
            coefficients[n] = sqrt(normalisedDoppler) * gam(0.75) / gam(1.25);

        coefficients[n] *= weight;
    }

    return coefficients;
}

// class OFDMChannel
OFDMChannel::OFDMChannel (const Model model, const double speed_kmh, const double frequency_Hz) : doppler(speed_kmh * frequency_Hz / (3.6 * 299792458.0)), _time(0.0)
{
    switch ( model )
    {
    case WINNER_SISO_URBAN_MACRO :
        _tapDelays .resize(18);
        _tapWeights.resize(18);

        for (size_t path = 0; path < 6; path++)
        {
            double pathWeight_dB;
            double pathDelay_s;
            switch ( path )
            {
            case 0:
                pathWeight_dB = 0.0;
                pathDelay_s   = 0.0;
                break;
            case 1:
                pathWeight_dB = -1.425;
                pathDelay_s   = 467.0E-9;
                break;
            case 2:
                pathWeight_dB = -4.217;
                pathDelay_s   = 1127.0E-9;
                break;
            case 3:
                pathWeight_dB = -7.852;
                pathDelay_s   = 1981.0E-9;
                break;
            case 4:
                pathWeight_dB = -12.037;
                pathDelay_s   = 3031.0E-9;
                break;
            case 5:
                pathWeight_dB = -14.919;
                pathDelay_s   = 4908.0E-9;
                break;
            }


            for (size_t subpath = 0; subpath < 3; subpath++)
            {
                double subpathWeight_dB;
                double subpathDelay_s;
                switch ( subpath )
                {
                case 0:
                    subpathWeight_dB = -3.031;
                    subpathDelay_s   = 0.0E-9;
                    break;
                case 1:
                    subpathWeight_dB = -5.229;
                    subpathDelay_s   = 7.0E-9;
                    break;
                case 2:
                    subpathWeight_dB = -6.990;
                    subpathDelay_s   = 27.0E-9;
                    break;
                }

                _tapWeights[path*3 + subpath] = pow(10.0, 0.05 * (pathWeight_dB + subpathWeight_dB));
                _tapDelays [path*3 + subpath] = pathDelay_s + subpathDelay_s;
            }
        }
        break;

    case WINNER_SISO_URBAN_MICRO:
        _tapDelays .resize(24);
        _tapWeights.resize(24);

        for (size_t path = 0; path < 6; path++)
        {
            double pathWeight_dB;
            double pathDelay_s;
            switch ( path )
            {
            case 0:
                pathWeight_dB = 0.0;
                pathDelay_s   = 0.0;
                break;
            case 1:
                pathWeight_dB = -0.783;
                pathDelay_s   = 261.0E-9;
                break;
            case 2:
                pathWeight_dB = -2.775;
                pathDelay_s   = 429.0E-9;
                break;
            case 3:
                pathWeight_dB = -4.605;
                pathDelay_s   = 608.0E-9;
                break;
            case 4:
                pathWeight_dB = -5.513;
                pathDelay_s   = 811.0E-9;
                break;
            case 5:
                pathWeight_dB = -7.658;
                pathDelay_s   = 1019.0E-9;
                break;
            }


            for (size_t subpath = 0; subpath < 4; subpath++)
            {
                double subpathWeight_dB;
                double subpathDelay_s;
                switch ( subpath )
                {
                case 0:
                    subpathWeight_dB = -4.559;
                    subpathDelay_s   = 0.0E-9;
                    break;
                case 1:
                    subpathWeight_dB = -6.021;
                    subpathDelay_s   = 5.0E-9;
                    break;
                case 2:
                    subpathWeight_dB = -6.990;
                    subpathDelay_s   = 11.0E-9;
                    break;
                case 3:
                    subpathWeight_dB = -6.990;
                    subpathDelay_s   = 28.0E-9;
                    break;
                }

                _tapWeights[path*4 + subpath] = pow(10.0, 0.05 * (pathWeight_dB + subpathWeight_dB));
                _tapDelays [path*4 + subpath] = pathDelay_s + subpathDelay_s;
            }
        }
        break;
    }

    // normalise tapWeights
    _tapWeights / _tapWeights.sum_of_squares();

    // initialise Rayleigh fading taps
    for (size_t i=0; i<_tapDelays.size(); i++)
        _taps.push_back(new RayleighFadingProcess(doppler, _tapWeights[i], _tapDelays[i]));
}

double OFDMChannel::frequency_response (const double frequency_Hz) const
{
    complex<double> voltageGain = complex<double>(0.0,0.0);
    for (size_t tap = 0; tap < _taps.size(); tap++)
        voltageGain += _taps[tap]->gain() * polar(1.0, -2.0 * M_PI * frequency_Hz * _taps[tap]->delay);

    return norm(voltageGain);
}

dVector OFDMChannel::frequency_response (const double startFrequency_Hz, const double stepFrequency_Hz, const size_t numFrequencies) const
{
    dVector powerGains(numFrequencies);

    for (size_t index = 0; index < numFrequencies; index++)
        powerGains[index] = frequency_response(startFrequency_Hz + stepFrequency_Hz * index);

    return powerGains;
}

dVector OFDMChannel::frequency_response (const double frequencyInterval_Hz, const double timeInterval_s) const
{
    size_t numPoints = round(1.0 / (frequencyInterval_Hz * timeInterval_s));

    // sample impulse response with resolution timeInterval_s
    cVector impulseResponse(numPoints, complex<double>(0.0,0.0));

    for (size_t tap = 0; tap < _taps.size(); tap++)
    {
        size_t bin = round(_taps[tap]->delay / timeInterval_s);
        if ( impulseResponse[bin] == complex<double>(0.0,0.0) )
            impulseResponse[bin] = _taps[tap]->gain();
        else
            impulseResponse[bin] = sqrt(sqr(impulseResponse[bin]) + sqr(_taps[tap]->gain()));
    }

    cVector fft = FFT::fft(impulseResponse, true);

    return norm(fft);
}

void OFDMChannel::set_time (const double time)
{
    _time = time;

    for (size_t i=0; i<_taps.size(); i++)
        _taps[i]->update(time);
}

OFDMChannel::RayleighFadingProcess::RayleighFadingProcess (const double doppler, const double weight, const double delay, const size_t numSinusoids)
: _numSinusoids(numSinusoids), _amplitude(weight / sqrt(static_cast<double>(_numSinusoids))), delay(delay), _time(-1.0)
{
    _frequencies.resize(numSinusoids, 0.0);
    for (size_t i=0; i<_numSinusoids; i++)
        _frequencies[i] = -2.0 * M_PI * doppler * cos(globalRandom.angle());

    _phaseOffsets.resize(numSinusoids, 0.0);
    for (size_t i=0; i<_numSinusoids; i++)
        _phaseOffsets[i] = globalRandom.angle();

    update(0.0);
}

void OFDMChannel::RayleighFadingProcess::update (const double time)
{
    if (time != _time)
    {
        _time = time;
        _gain = complex<double>(0.0,0.0);
        for (size_t i=0; i<_numSinusoids; i++)
        {
            const double phase = _frequencies[i]*time - _phaseOffsets[i];
            _gain += polar(_amplitude, phase);
        }
    }
}

/*
cVector MultipathRayleighChannel::operator()(const cVector &input, cMatrix &channel)
{
    const size_t length   = input.size();
    const size_t maxDelay = tapDelays.back();
    cVector output(length + maxDelay);

    const size_t fftSize         = pow(2.0, ceil(log2(length)));
    const size_t numNoiseSamples = ceil(normalisedDoppler * fftSize);

    dVector F(fftSize);
    F.ramp(0.0, 1.0 / fftSize, 1.0, 0.5);

    cVector S(fftSize, 0.0);
    double  norm = 0.0;
    for (size_t i=0; i<fftSize; i++)
        if ( fabs(F[i]) < normalisedDoppler )
            {
                const double temp = 1.5 / (M_PI * normalisedDoppler * sqrt(1.0-sqr(F[i]/normalisedDoppler)));
                norm += temp;
	            S[i] = sqrt(temp) * fftSize;
            }
	
    S /= norm;
	
    channel.set_size(tapDelays.size(), length, 0.0);
    for (size_t i=0; i<tapDelays.size(); i++)
    {
        cVector randomSamples;
        for (size_t j=0; j<numNoiseSamples; j++)
            randomSamples.push_back(globalRandom.complex_gaussian());
        for (size_t j=0; j<fftSize - 2*numNoiseSamples; j++)
            randomSamples.push_back(0.0);
        for (size_t j=0; j<numNoiseSamples; j++)
            randomSamples.push_back(globalRandom.complex_gaussian());

        cVector x = FFT::ifft(S * randomSamples);
        cVector y = x.mid(0,length);

        for (size_t j=0; j<length; j++)
            channel[i][j] = y[j] * normalisedLinearTapWeigths[i];

        cVector tapContribution(length + maxDelay);
        tapContribution.replace_mid(tapDelays[i], input * channel[i]);

        output += tapContribution;
    }
	
    return AWGNChannel::operator ()(output);
}
*/