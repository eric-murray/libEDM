#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_types.h>
#include <libEDM_library.h>
#include <libEDM_random.h>

class ShadowFading {
public:
    const dB standardDeviation;
    const dB median; 
    const dB maxDeviation;

   	operator double() const {return _shadowFading;}

    ShadowFading (const dB &standardDeviation, const dB median = dB(0.0), const dB maxDeviation = dB(dINFINITY)) : standardDeviation(standardDeviation), median(median), maxDeviation(abs(maxDeviation)), _shadowFading(generateFadeValue()) {assert( this->maxDeviation >= 0.0 );}

    void reset() {if ( standardDeviation != 0.0 ) _shadowFading = generateFadeValue();}

private:
    double _shadowFading; // linear fade value

    double generateFadeValue() const;
};

class Propagation {
public:
	enum ClutterType {Default, DenseUrban, Urban, SubUrban, Rural};

    const MHz         frequency;		// MHz
    const metres      minDistance;		// m
	const ClutterType defaultClutterType;

	Propagation (const MHz &frequency_MHz, const ClutterType defaultClutterType, const metres &minDistance = metres(0.0))
		: frequency         (frequency_MHz), 
          _logFrequency     (log10(frequency_MHz)),
          _logFrequencyLimit(std::min(std::max(log10(150.0), _logFrequency), log10(2000.0))),
          minDistance       (minDistance),
          _defaultk         (-1.0),
          _1_defaultk       (-1.0),
          _defaultGamma     (-1.0),
          defaultClutterType(defaultClutterType)
    {}

	Propagation (const MHz &frequency_MHz, const dB &defaultkdB, const double defaultGamma, const ClutterType defaultClutterType = Urban, const metres &minDistance = metres(0.0))
		: frequency         (frequency_MHz),
          _logFrequency     (log10(frequency_MHz)),
          _logFrequencyLimit(std::min(std::max(log10(150.0), _logFrequency), log10(2000.0))),
          minDistance       (minDistance),
          _defaultk         (dB2linear(defaultkdB)),
          _1_defaultk       (1.0/_defaultk),
          _defaultGamma     (defaultGamma),
          defaultClutterType(defaultClutterType)
    {}

    double pathgain_LOS          (const metres &distance_m) const;
    double pathgain_power_law    (const metres &distance_m) const;
    double pathgain_power_law    (const metres &distance_m, const dB &k_dB, const double gamma) const;
	double pathgain_ExtendedHata (const metres &distance_m, const metres &bsHeight_m, const metres &msHeight_m, const ClutterType clutterType, bool &isLOSLimited) const;
    double pathgain_ExtendedHata (const metres &distance_m, const metres &bsHeight_m, const metres &msHeight_m, const ClutterType clutterType = Default) const;

private:
    const static double _speedOfLight;

	const double _defaultk, _1_defaultk;
	const double _defaultGamma;
	const double _logFrequency, _logFrequencyLimit;

	double pathloss_ExtendedHata_dB (const metres &distance_m, const metres &bsHeight_m, const metres &msHeight_m, const ClutterType clutterType, const bool applyLOSLimit = true) const;

    double HataA (const metres &msHeight_m, const ClutterType clutterType) const;
	double HataB (const metres &bsHeight_m) const;
};