#define _USE_MATH_DEFINES
#include <complex>
#include <cmath>
#include <string>

#include <libEDM_fft.h>
#include <libEDM_library.h>
#include <libEDM_random.h>

using std::abs;
using std::conj;
using std::min;
using std::max;
using std::norm;
using std::pair;
using std::string;

void error (string &errorMessage)
{
	cerr << errorMessage << endl;
	abort();
}

void error (char *errorMessage)
{
	error(string(errorMessage));
}

cVector conj (const cVector &x)
{
    cVector output;
    for (size_t i=0; i<x.size(); i++)
        output.push_back(conj(x[i]));
    return output;
}

dVector norm (const cVector &x)
{
    dVector output;
    for (size_t i=0; i<x.size(); i++)
        output.push_back(norm(x[i]));
    return output;
}

dVector real (const cVector &x)
{
    dVector output;
    for (size_t i=0; i<x.size(); i++)
        output.push_back(x[i].real());
    return output;
}

dVector imag (const cVector &x)
{
    dVector output;
    for (size_t i=0; i<x.size(); i++)
        output.push_back(x[i].imag());
    return output;
}

dVector abs (const cVector &x)
{
    dVector output;
    for (size_t i=0; i<x.size(); i++)
        output.push_back(abs(x[i]));
    return output;
}

dVector phase (const cVector &x)
{
    dVector output;
    for (size_t i=0; i<x.size(); i++)
        output.push_back(atan2(x[i].real(), x[i].imag()));
    return output;
}

dVector xcorr (const dVector &x, const dVector &y)
{
    const size_t vectorLength = x.size();

    // copy dVectors to cVectors
    cVector cx(vectorLength);
    cVector cy(vectorLength);
    for (size_t i=0; i<vectorLength; i++)
    {
        cx[i] = x[i];
        cy[i] = y[i];
    }

    cVector ccorr = xcorr(cx, cy);

    // copy real parts to output
    dVector output;
    for (size_t i=0; i<2*vectorLength; i++)
        output.push_back(ccorr[i].real());

    return output;
}

cVector xcorr (const cVector &x, const cVector &y)
{
    const size_t vectorLength = x.size();

    // check vectors are same length
    assert( vectorLength == y.size() );

    // zero pad vectors
    cVector cx = x;
    cVector cy = y;
    cx.resize(2*vectorLength);
    cy.resize(2*vectorLength);

    // compute ffts
    cVector fftx = FFT::fft(cx);
    cVector ffty = FFT::fft(cy);

    // take inverse of multiplication to get xcorrelation
    cVector ccorr = conj(FFT::ifft(fftx * conj(ffty)));

    // re-sort output, most negative lag first
    cVector output = ccorr.right(vectorLength);
    output.ins(vectorLength, ccorr.left(vectorLength));

    return output;
}

size_t combinations (const size_t N, const size_t k)
{
    if ( (N==k) || (k==0) )
        return 1;
    
    size_t answer = max(N-k+1, k+1);
    size_t num  = answer+1;
    size_t den  = 2;

    while (num <= N)
    {
        answer *= num++;
        answer /= den++;
    }

    return answer;
}

double binomial (const size_t N, const size_t k, const double p)
{
    return combinations(N,k) * pow(p,k) * pow((1.0-p),(N-k));
}

double binomial_CDF (const size_t N, const size_t y, const double p)
{
    double answer = 0.0;
    for (size_t k=0; k<=y; k++)
        answer += binomial(N,k,p);

    return answer;
}

//
// class LookUpTable
//

LookUpTable::LookUpTable (ifstream &inputFile, const double yMinLimit, const double yMaxLimit) : _yMinLimit(yMinLimit), _yMaxLimit(yMaxLimit), _minX(dINFINITY), _minY(dINFINITY), _maxX(-dINFINITY), _maxY(-dINFINITY)
{
    if ( !inputFile )
        // NULL input file passed
        error("NULL input file passed to LookUpTable constructor");

    while ( !inputFile.eof() )
    {
        // read x value
        char xBuffer[30];
        inputFile.get(&xBuffer[0], 30, ',');

        // read comma delimiter
        inputFile.get();

        // read y value
        char yBuffer[30];
        inputFile.get(&yBuffer[0], 30);

        // read newline character
        inputFile.get();

        string xString = xBuffer;
        string yString = yBuffer;

        if ( !(xString.empty() || yString.empty()) )
            add_point(atof(&xBuffer[0]), atof(&yBuffer[0]));
    }
}

void LookUpTable::add_point (const double x, double y)
{
    y = std::max(y, _yMinLimit);
    y = std::min(y, _yMaxLimit);

    _points.insert(pair<double,double>(x,y));

    _minX = min(_minX, x);
    _minY = min(_minY, y);
    _maxX = max(_maxX, x);
    _maxY = max(_maxY, y);

    // check if LUT points are regularly sampled
    _regularXSampling = true;
    _samplingInterval = -1.0;

    double currentX = 0.0;
    for (const_iterator point = _points.begin(); point != _points.end(); point++)
    {
        double previousX = currentX;
        currentX = point->first;

        if (point != _points.begin())
        {
            double samplingInterval = currentX - previousX;

            if (_samplingInterval == -1.0)
                _samplingInterval = samplingInterval;
            else
                if (_samplingInterval != samplingInterval)
                {
                    _regularXSampling = false;
                    break;
                }
        }
    }

    // store y values in vector if regularly sampled
    if ( _regularXSampling )
    {
        _pointVector.clear();
        for (const_iterator point = _points.begin(); point != _points.end(); point++)
            _pointVector.push_back(point->second);
    }

    // check if monotonic rising or falling
    _monotonicRising  = true;
    _monotonicFalling = true;

    double currentY = 0.0;
    for (const_iterator point = _points.begin(); point != _points.end(); point++)
    {
        double previousY = currentY;
        currentY = point->second;

        if (point != _points.begin())
        {
            if (currentY > previousY)
                _monotonicFalling = false;

            if (currentY < previousY)
                _monotonicRising = false;
        }
    }
}

void LookUpTable::print()
{
    for (const_iterator point = _points.begin(); point != _points.end(); point++)
    {
        cout << point->first << ",";
        cout << point->second << endl;
    }
    cout << endl;
}

double LookUpTable::cdf_sample () const
{
    return x(globalRandom.sample());
}

double LookUpTable::random_x () const
{
    return x(globalRandom.sample() * (_maxY - _minY) + _minY);
}

double LookUpTable::random_y () const
{
    return y(globalRandom.sample() * (_maxX - _minX) + _minX);
}

double LookUpTable::x (const double y) const
{
    const_iterator upper = _points.begin(), lower;
    for (const_iterator point = _points.begin(); point != _points.end(); point++)
    {
        lower = upper;
        upper = point;

        if ( upper != _points.begin() )
        {
            if ( _monotonicFalling && (upper->second < y) )
                break;
            if ( !_monotonicFalling && (upper->second > y) )
                break;
        }
    }

    double xDiff = upper->first  - lower->first;
    double yDiff = upper->second - lower->second;
    double x = lower->first + xDiff * ((y - lower->second) / yDiff);

    return x;
}

int LookUpTable::x (const int y) const
{
    return static_cast<int>(x(static_cast<double>(y)));
}

double LookUpTable::y (const double x) const
{
    double y;
    if ( _regularXSampling )
    {
        // LUT is regularly sampled, so use _pointVector
        int pointId = static_cast<int>((x - _points.begin()->first) / _samplingInterval);

        pointId = max(pointId, 0);
        pointId = min(pointId, static_cast<int>(numPoints()) - 2);

        double upper = _pointVector[pointId + 1];
        double lower = _pointVector[pointId];
        y = lower + (upper - lower) * (x - (pointId * _samplingInterval + _points.begin()->first)) / _samplingInterval;
    }
    else
    {
        // LUT is not regularly sampled
        const_iterator upper = _points.begin(), lower;
        for (const_iterator point = _points.begin(); point != _points.end(); point++)
        {
            lower = upper;
            upper = point;

            if ( (upper != _points.begin()) && (upper->first > x) )
                break;
        }

        double xDiff = upper->first  - lower->first;
        double yDiff = upper->second - lower->second;
        y = lower->second + yDiff * ((x - lower->first) / xDiff);
    }

    y = std::max(y, _yMinLimit);
    y = std::min(y, _yMaxLimit);

    return y;
}

int LookUpTable::y(const int x) const
{
    return static_cast<int>(y(static_cast<double>(x)));
}