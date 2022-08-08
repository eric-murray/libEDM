#define _USE_MATH_DEFINES
#include <complex>
#include <math.h>

#include <libEDM_filter.h>
#include <libEDM_library.h>

dVector RRCFilter::impulseResponse (const size_t filterLength, const double alpha, const size_t overSamplingFactor)
{
    dVector impulseResponse(filterLength * overSamplingFactor + 1);

    for (size_t i=0; i<impulseResponse.size(); i++)
    {
        const double t = (i - 0.5 * overSamplingFactor * filterLength) / overSamplingFactor;

        // 3GPP pulse shaping filter
        double tapWeight = ( sin(M_PI*t*(1.0-alpha)) + (4.0*alpha*t)*cos(M_PI*t*(1.0+alpha)) ) / ( M_PI*t*(1.0-sqr(4.0*alpha*t)) );
        impulseResponse[i] = tapWeight;
    }

    // fix for t=0
    impulseResponse[overSamplingFactor*filterLength/2] = 4.0*alpha/M_PI + 1.0 - alpha;

    // need fix for t = +/- 1.0/(4.0 * alpha)
    // this only occurs if overSamplingFactor is an integer multiple of 100*alpha (unlikely)

    // scale by overSamplingFactor
    impulseResponse /= sqrt(static_cast<double>(overSamplingFactor));

    return impulseResponse;
}

