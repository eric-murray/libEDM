#include "libEDM_histogram.h"

using std::ios_base;
using std::setiosflags;

//
// class LogHistogram
//

string LogHistogram::dBType() const
{
    switch (_logType)
    {
	case B:
		return string("Bells");
    case dB:
        return string("dB");
    case dBm:
        return string("dBm");
    case dBW:
        return string("dBW");
    default:
		error("LogHistogram::dBType logType not supported");
    }
}

void LogHistogram::print_footer(ostream &file) const
{
	file << setiosflags(ios_base::scientific);

	file << "Mean Value (dB)";
    for (size_t i = 0; i < 4; i++)
        file << ", " << linear2dB(mean());
    file << endl;

    file << "Std. Dev.  (dB)";
    for (size_t i = 0; i < 4; i++)
        file << ", " << linear2dB(stdDev());
    file << endl << endl;

    file << "Logarithmic Mean Value (" << dBType() << ")";
    for (size_t i = 0; i < 4; i++)
        file << ", " << dBMean();
    file << endl;

    file << "Logarithmic Std. Dev.  (" << dBType() << ")";
    for (size_t i = 0; i < 4; i++)
        file << ", " << dBStdDev();
    file << endl << endl;

    file.flush();
}

void LogHistogram::reset()
{
    Histogram<double>::reset();

    _dBSum          = 0.0;
    _dBSumOfSquares = 0.0;
}

void LogHistogram::update (const double sample)
{
    Histogram<double>::update(sample);

    double logSample;
    switch (_logType)
    {
	case B:
		logSample = log10(sample);
		break;
	
	case dB:
    case dBW:
        logSample = linear2dB(sample);
        break;

    case dBm:
        logSample = watts2dBm(sample);
        break;
    }

    _dBSum          += logSample;
    _dBSumOfSquares += sqr(logSample);

	update_histogram(logSample);
}