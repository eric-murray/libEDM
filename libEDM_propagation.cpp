#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_library.h>
#include <libEDM_propagation.h>

using std::max;
using std::min;
using std::pow;

//
// class ShadowFading
//
double ShadowFading::generateFadeValue () const
{
    dB shadowFading = dB(globalRandom.gaussian(median, standardDeviation));

    shadowFading = std::min(shadowFading, median + maxDeviation);
    shadowFading = std::max(shadowFading, median - maxDeviation);

    return dB2linear(shadowFading);
}

//
// class Propagation
//

const double Propagation::_speedOfLight = 299792458.0; // m/s

double Propagation::HataA (const metres &msHeight_m, const ClutterType clutterType) const
{
	switch (clutterType)
	{
	case DenseUrban:
		if ( frequency < MHz(300.0) )
			return 8.29 * sqr(log10(msHeight_m * 1.54)) - 1.1;
		else
			return 3.2 * sqr(log10(msHeight_m * 11.75)) - 4.97;

	case Urban:
		return (1.1 * _logFrequency - 0.7) * msHeight_m - (1.56 * _logFrequency - 0.8);

	default:
		error("Clutter type not defined");

	}
}

double Propagation::HataB (const metres &bsHeight_m) const
{
	// 20.0 * log10(30.0) = 29.542425094393248745900558065102
	return min(0.0, 20.0*log10(bsHeight_m) - 29.542425094393248745900558065102);
}

double Propagation::pathgain_ExtendedHata (const metres &distance_m, const metres &bsHeight_m, const metres &msHeight_m, const ClutterType clutterType) const
{
    bool dummy;

    return pathgain_ExtendedHata(distance_m, bsHeight_m, msHeight_m, clutterType, dummy);
}

double Propagation::pathgain_ExtendedHata (const metres &distance_m, const metres &bsHeight_m, const metres &msHeight_m, const ClutterType clutterType, bool &isLOSLimited) const
{
	const metres distance = max(distance_m, minDistance);

    double pathgain = pow(10.0, -0.1*pathloss_ExtendedHata_dB(distance, bsHeight_m, msHeight_m, clutterType == Default ? defaultClutterType : clutterType, false));

    const double losPathgain = pathgain_LOS(distance_m);
    if ( pathgain > losPathgain )
    {
        pathgain = losPathgain;
        isLOSLimited = true;
    }
    else
        isLOSLimited = false;

    return pathgain;
}

double Propagation::pathgain_LOS (const metres &distance_m) const
{
	const metres distance = max(distance_m, minDistance);

    return sqr(_speedOfLight / (4.0 * M_PI * distance * frequency * 1.0E6));
}

double Propagation::pathgain_power_law (const metres &distance_m) const
{
	const metres distance = max(distance_m, minDistance);

	return _1_defaultk * pow(distance, -_defaultGamma);
}

double Propagation::pathgain_power_law (const metres &distance_m, const dB &k_dB, const double gamma) const
{
	const metres distance = max(distance_m, minDistance);

    return dB2linear(-k_dB) * pow(distance, -gamma);
}

double Propagation::pathloss_ExtendedHata_dB (const metres &distance_m, const metres &bsHeight_m, const metres &msHeight_m, const ClutterType clutterType, const bool applyLOSLimit) const
{
	if ( frequency < MHz(30.0) )
		error("Extended Hata model is not valid below 30 MHz");

	if ( frequency > MHz(3000.0) )
		error("Extended Hata model is not valid above 3000 MHz");

	const metres maxBSHeight_m = metres(30.0);

    double pathloss_dB;
	switch (clutterType)
	{
	case DenseUrban:
	case Urban:
		{
			const double logBSHeight  = log10(max(bsHeight_m, maxBSHeight_m));
			const double distanceTerm = (44.9 - 6.55*logBSHeight) * log10(distance_m * 1.0E-3);

			if ( frequency <= 150.0 )
				// 26.16 * log10(150.0) = 56.93
				pathloss_dB = 69.55 + 26.16*_logFrequency - 13.82*logBSHeight - HataA(msHeight_m, clutterType) - HataB(bsHeight_m) + distanceTerm + 56.93 - 20.0*log10(150.0 / frequency);
			
			else if ( frequency <= 1500.0 )
				pathloss_dB = 69.55 + 26.16*_logFrequency - 13.82*logBSHeight - HataA(msHeight_m, clutterType) - HataB(bsHeight_m) + distanceTerm;

			else if ( frequency <= 2000.0 )
				pathloss_dB = 46.3 + 33.9*_logFrequency - 13.82*logBSHeight - HataA(msHeight_m, Urban) - HataB(bsHeight_m) + distanceTerm;

			else if ( frequency <= 3000.0 )
				// 33.9 * log10(2000) - 10 * log10(2000) = 78.89
				pathloss_dB = 46.3 + 78.89 + 10.0*_logFrequency - 13.82*logBSHeight - HataA(msHeight_m, Urban) - HataB(bsHeight_m) + distanceTerm;

            else
        		error("Coding error - this point should not be reached");
		}
		break;

	case SubUrban:
	case Rural:
		{
			switch (clutterType)
			{
			case SubUrban:
                // log10(28.0) = 1.447
				pathloss_dB = pathloss_ExtendedHata_dB(distance_m, bsHeight_m, msHeight_m, Urban, false) - 2.0 * sqr(_logFrequencyLimit - 1.447) - 5.4;
                break;

			case Rural:
				pathloss_dB = pathloss_ExtendedHata_dB(distance_m, bsHeight_m, msHeight_m, Urban, false) - 4.78 * sqr(_logFrequencyLimit) + 18.33 * _logFrequencyLimit - 40.94;
                break;

            default:
        		error("Coding error - this point should not be reached");
			}
		}
        break;

	default:
		error("Coding error - this point should not be reached");
	}

    if ( applyLOSLimit )
    {
        double losPathloss_dB = -27.6 + 20.0 * log10(distance_m) + 20.0 * _logFrequency;
        
        if ( pathloss_dB < losPathloss_dB )
            pathloss_dB = losPathloss_dB;
    }

    return pathloss_dB;
}
