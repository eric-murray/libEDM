#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "libEDM_types.h"
#include "libEDM_matrix.h"

namespace fs = boost::filesystem;

using std::ceil;
using std::cerr;
using std::cout;
using std::endl;
using std::fixed;
using std::ifstream;
using std::log10;
using std::max;
using std::multimap;
using std::numeric_limits;
using std::ofstream;
using std::ostream;
using std::pow;
using std::scientific;
using std::setprecision;
using std::setw;
using std::sqrt;
using std::string;
using std::vector;

const double dINFINITY   = numeric_limits<double>::infinity();
const size_t uMAX        = numeric_limits<size_t>::max();
const double M_1_LOG10_2 = 3.3219280948873623478703194294894;

__declspec(noreturn) void error(string &errorMessage);
__declspec(noreturn) void error(char   *errorMessage);

cVector xcorr (const cVector &x, const cVector &y);
dVector xcorr (const dVector &x, const dVector &y);
cVector conj  (const cVector &x);
dVector norm  (const cVector &x);
dVector real  (const cVector &x);
dVector imag  (const cVector &x);
dVector abs   (const cVector &x);
dVector phase (const cVector &x);

size_t combinations (const size_t N, const size_t k);
double binomial     (const size_t N, const size_t k, const double p);
double binomial_CDF (const size_t N, const size_t y, const double p);

inline double pow   (const double x, const size_t y) {return pow(x, static_cast<int>(y));}
inline int    round (const double x)                 {return static_cast<int>(floor(x + 0.5));}
inline double sgn   (const double x)                 {return x==0.0 ? 0.0 : (x<0.0 ? -1.0 : 1.0);}

// templated trigonometric functions
template <class Type> inline Type atan (const double x)                 {return Type(atan(x));}
template <class Type> inline Type atan2(const double x, const double y) {return Type(atan2(x,y));}
template <class Type> inline Type acos (const double x)                 {return Type(acos(x));}
template <class Type> inline Type asin (const double x)                 {return Type(asin(x));}

template <class Type> inline Type   abs           (const Type &x)                {return x < 0.0 ? -x : x;}
template <class Type> inline Type   sqr           (const Type &x)                {return x*x;}
template <class Type> inline Type   sqrt          (const Type &x)                {return Type(sqrt(x()));}
template <class Type> inline double divide        (const Type &x, const Type &y) {return static_cast<double>(x) / static_cast<double>(y);}
template <class Type> inline Type   most_negative ()                             {return numeric_limits<Type>::is_integer ? numeric_limits<Type>::min() : -numeric_limits<Type>::infinity();} 

//
// Value conversion routines
//

inline double dB2linear (const dB  x) {return pow(10.0, 0.1*x);}
inline double dBm2Watts (const dBm x) {return dB2linear(x - dB(30.0));}

inline dB linear2dB (double linear, double minValue = -numeric_limits<double>::infinity())
{
    if (linear > 0.0)
        return dB(max(minValue, 10.0*log10(linear)));
    else
        return dB(minValue);
}

inline dBm watts2dBm (double watts, double minValue = -numeric_limits<double>::infinity())
{
    return dBm(linear2dB(watts, minValue-30.0) + 30.0);
}

//
// class LookUpTable
//

class LookUpTable {
public:
    size_t numPoints() const {return static_cast<size_t>(_points.size());}

    LookUpTable (                     const double yMinLimit = -dINFINITY, const double yMaxLimit = dINFINITY) : _yMinLimit(yMinLimit), _yMaxLimit(yMaxLimit), _minX(dINFINITY), _minY(dINFINITY), _maxX(-dINFINITY), _maxY(-dINFINITY) {}
    LookUpTable (ifstream &inputFile, const double yMinLimit = -dINFINITY, const double yMaxLimit = dINFINITY);

    double x (const double y) const;
    double y (const double x) const;
    int    x (const int    y) const;
    int    y (const int    x) const;

    void   add_point (const double x, double y);
    double random_x  () const;
    double random_y  () const;
    double cdf_sample() const;

    void print();

private:
    typedef multimap<double,double>::const_iterator const_iterator;

    const double _yMinLimit, _yMaxLimit;

    bool                    _regularXSampling, _monotonicRising, _monotonicFalling;
    double                  _samplingInterval;
    double                  _minX, _minY, _maxX, _maxY;
    multimap<double,double> _points;
    vector<double>          _pointVector;
};