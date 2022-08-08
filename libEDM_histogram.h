#pragma once

#include "libEDM_library.h"

//
// class Histogram
//
struct DataType
{
    long   bin;
    size_t numSamples;
    double pdf;
    double cdf;

    double inverseCDF() {return 1.0 - cdf;}
};

typedef vector<DataType> HistogramData;

template <class Type>
class Histogram {
public:
    const string title;
    const string columnLabel;
    const string defaultFilename;

    size_t count()        const {return _count;}
    Type   sum()          const {return _sum;}
    Type   sumOfSquares() const {return _sumOfSquares;}

    Type   high()         const {return _low + _histogram.size() * interval;}
	Type   interval()     const {return _interval;}
	Type   max()          const {return _max;}
    Type   min()          const {return _min;}
    double mean()         const {return _count == 0 ? 0.0 : _sum / _count;}
    double meanSquare()   const {return _count == 0 ? 0.0 : _sumOfSquares / _count;}
    double stdDev()       const {return sqrt(variance());}
    double variance()     const {return meanSquare() - sqr(mean());}

    // constructors
    Histogram (const string &title, const string &columnLabel, Type interval, size_t maxBins, const string &defaultFilename)
        : title(title), columnLabel(columnLabel), _defaultInterval(interval), _maxBins(maxBins), defaultFilename(defaultFilename) {reset();}

    Type getPercentile (double percentile) const;

    void reset ();

    // print methods
    void print (ostream        &file,            const Type startFrom = most_negative<Type>(), const Type stopAt = numeric_limits<Type>::infinity(), const bool suppressEmptyBins = false) const;
	void print (const char     *fileName = NULL, const Type startFrom = most_negative<Type>(), const Type stopAt = numeric_limits<Type>::infinity(), const bool suppressEmptyBins = false) const;
	void print (const string   &fileName,        const Type startFrom = most_negative<Type>(), const Type stopAt = numeric_limits<Type>::infinity(), const bool suppressEmptyBins = false) const {print(fileName.c_str(), startFrom, stopAt, suppressEmptyBins);}
    void print (const fs::path &file,            const Type startFrom = most_negative<Type>(), const Type stopAt = numeric_limits<Type>::infinity(), const bool suppressEmptyBins = false) const {print(file.string(),    startFrom, stopAt, suppressEmptyBins);}
    void print (                                 const Type startFrom                        , const Type stopAt = numeric_limits<Type>::infinity(), const bool suppressEmptyBins = false) const {print(NULL            , startFrom, stopAt, suppressEmptyBins);}

    // pure virtual functions
    // these stop the Histogram class being instantiate directly
    virtual void update(const Type sample) = 0;

protected:
	const Type _defaultInterval;

    size_t         _count, _maxBins, _positiveInfinities, _negativeInfinities;
    Type           _interval, _low, _max, _min;
	double         _sum, _sumOfSquares;
    vector<size_t> _histogram;

    long which_bin (const Type sample) const;

    Type binCentre (const long bin) const;

	HistogramData get_data (const bool suppressEmptyBins = false) const;

    void update_histogram (const Type sample);

	virtual void print_footer(ostream &file) const {};
};

template <class Type>
void Histogram<Type>::update(const Type sample)
{
    _count++;

    if (sample < _min)
        _min = sample;

    if (sample > _max)
        _max = sample;

    _sum          += static_cast<double>(sample);
    _sumOfSquares += sqr(static_cast<double>(sample));
}

template <class Type>
HistogramData Histogram<Type>::get_data (const bool suppressEmptyBins) const
{
    HistogramData data;

	if ( _negativeInfinities > 0 )
	{
		const double value = divide(_negativeInfinities, _count);

        DataType binData;
        binData.bin        = numeric_limits<long>::min();
        binData.numSamples = _negativeInfinities;
        binData.pdf        = value;
        binData.cdf        = value;

        data.push_back(binData);
	}

    size_t cumulative = _negativeInfinities;
    long   bin        = 0;

	// skip bins with no hits
	while ( _histogram[bin] == 0 )
		bin++;

	if ( !suppressEmptyBins && _negativeInfinities == 0 && !(numeric_limits<Type>::is_modulo && binCentre(bin) == numeric_limits<Type>::min()) )
	{
		// store bin before first bin to have hits unless Type is wraparound and bin centre of lowest bin with hits is equal to minimum value
        DataType binData;
        binData.bin        = bin - 1;
        binData.numSamples = 0;
        binData.pdf        = 0.0;
        binData.cdf        = 0.0;

        data.push_back(binData);
	}

    do
    {
        cumulative += _histogram[bin];

		if ( !suppressEmptyBins || _histogram[bin] > 0 )
		{
            DataType binData;
            binData.bin        = bin;
            binData.numSamples = _histogram[bin];
            binData.pdf        = divide(_histogram[bin], _count);
            binData.cdf        = divide(cumulative,      _count);

            data.push_back(binData);
		}

		bin++;

    } while (cumulative + _positiveInfinities < _count);

	if ( _positiveInfinities > 0 )
	{
        DataType binData;
        binData.bin        = numeric_limits<long>::max();
        binData.numSamples = _positiveInfinities;
        binData.pdf        = divide(_positiveInfinities, _count);
        binData.cdf        = 1.0;

        data.push_back(binData);
	}

	assert ( cumulative + _positiveInfinities == _count );

	return data;
}

template <class Type>
Type Histogram<Type>::getPercentile (double percentile) const
{
    percentile *= 0.01;
    percentile = std::min(percentile, 1.0);
    percentile = std::max(percentile, 0.0);

    HistogramData data = this->get_data(true);

    for (size_t i=0; i < data.size(); i++)
        if ( data[i].cdf > percentile )
            return binCentre(data[i].bin);

    // should only get to this point if percentile = 100.0
    assert ( percentile == 1.0 );

    // percentile not found, so return last bin value
    return binCentre(data.back().bin);
}

template <class Type>
void Histogram<Type>::print(ostream &file, const Type startFrom, const Type stopAt, const bool suppressEmptyBins) const
{
    const size_t numColumns = 5;

	if ( ! file )
		// file not open
		error("File not open");

    file << title.c_str() << endl;
	if ( _count == 0 )
	{
		file << "No data" << endl;
		return;
	}

    HistogramData data = this->get_data(suppressEmptyBins);

    file << "Number of samples";
    for (size_t i=1; i<numColumns; i++)
        file << ", " << _count;
    file << endl;

    file << "Mean Value";
    for (size_t i=1; i<numColumns; i++)
        file << ", " << scientific << mean();
    file << endl;

    file << "Std. Dev.";
    for (size_t i=1; i<numColumns; i++)
        file << ",  " << scientific << stdDev();
    file << endl;
    
    file << endl;

    this->print_footer(file);

    // print colums labels
    file << ", Num Samples, PDF, CDF, 1 - CDF" << endl;
    for (size_t i = 1; i < numColumns; i++)
        file << ", " << columnLabel;
    file << endl;

    const int startBin = which_bin(startFrom);
    const int stopBin  = which_bin(stopAt);

    // print entries from startFrom value up to first data entry
    if ( startFrom > most_negative<Type>() )
    {
        for (int bin = startBin; bin < data[0].bin; bin++)
        {
			file << setw(8) << fixed << setprecision(2) << binCentre(bin) << ",";
            file << setw(8) << fixed << 0 << ", ";
		    file << setw(9) << setprecision(4) << 0.0 << ",";
		    file << setw(9) << setprecision(4) << 0.0 << ",";
		    file << setw(9) << setprecision(4) << 1.0 << endl;
        }
    }

	for (size_t i=0; i<data.size(); i++)
	{
        // don't print values where x is less than startFrom value
        if ( data[i].bin < startBin || data[i].bin > stopBin )
            continue;

		if ( data[i].bin == numeric_limits<long>::min() )
			file << "  -1.#INF,";

		else if ( data[i].bin == numeric_limits<long>::max() )
			file << "   1.#INF,";

		else
			file << setw(8) << fixed << setprecision(2) << binCentre(data[i].bin) << ",";

		file << setw(9) << fixed << data[i].numSamples << ",";
		file << setw(9) << setprecision(4) << data[i].pdf << ",";
		file << setw(9) << setprecision(4) << data[i].cdf << ",";
		file << setw(9) << setprecision(4) << data[i].inverseCDF() << endl;
	}

    file << endl;

	file.flush();
}

template <class Type>
void Histogram<Type>::print (const char *inputFilename, Type startFrom, Type stopAt, const bool suppressEmptyBins) const
{
    const char *filename = inputFilename == NULL ? defaultFilename.c_str() : inputFilename;
	ofstream file(filename);

	if ( ! file )
		// Cannot open file
		error(string("Cannot open file ") + filename);

	print(file, startFrom, stopAt, suppressEmptyBins);
	file.close();
}


template <class Type>
void Histogram<Type>::reset()
{
    _count        = 0;
	_interval     = _defaultInterval;
	_low          = 0;
    _max          = most_negative<Type>();
    _min          = numeric_limits<Type>::max();
    _sum          = 0.0;
    _sumOfSquares = 0.0;

    _histogram.clear();
	_positiveInfinities = 0;
	_negativeInfinities = 0;

    // check if file with defaultFilename can be used
    if ( defaultFilename != "" )
    {
        ofstream file(defaultFilename.c_str());
	    if ( ! file )
		    // Cannot open file
		    error(string("Cannot open file ") + defaultFilename.c_str());
        else
            file.close();
    }
}

template <class Type>
Type Histogram<Type>::binCentre (const long bin) const
{
    if ( bin == numeric_limits<long>::min() )
        return most_negative<Type>();

    if ( bin == numeric_limits<long>::max() )
        return numeric_limits<Type>::infinity();

    return _low + bin * _interval;
}


template <class Type>
long Histogram<Type>::which_bin (const Type sample) const
{
	if ( numeric_limits<Type>::has_infinity )
	{
		// Test whether sample is +ve or -ve infinity
		if ( sample == numeric_limits<Type>::infinity() )
            return numeric_limits<long>::max();
		if ( sample == -numeric_limits<Type>::infinity() )
            return numeric_limits<long>::min();
	}
    return static_cast<long>( floor( (static_cast<double>(sample) - _low) / _interval ));
}

template <class Type>
void Histogram<Type>::update_histogram(const Type sample)
{
	if ( numeric_limits<Type>::has_infinity )
	{
		// Test whether sample is +ve or -ve infinity
		// Don't include these values directly in the histogram
		if ( sample == numeric_limits<Type>::infinity() )
		{
			_positiveInfinities++;
			return;
		}
		if ( sample == -numeric_limits<Type>::infinity() )
		{
			_negativeInfinities++;
			return;
		}
	}

	if ( _count == 1 )
		// set _low equal to quantised sample value
		_low = floor(divide(sample, _interval)) * _interval;
		//_low = ceil(divide(sample, _interval)) * _interval;

	long bin = which_bin(sample);
    if (bin < 0)
	{
		// insert bins at lower end of histogram
        _histogram.insert(_histogram.begin(), abs(bin), 0);

		// recalibrate _low
		_low = floor(divide(sample, _interval)) * _interval;
		//_low = ceil(divide(sample, _interval)) * _interval;

		bin = 0;
	}
    else
		if ( static_cast<size_t>(bin) >= _histogram.size() )
			_histogram.resize(bin+1, 0);

    _histogram[bin]++;

	if ( _histogram.size() > _maxBins )
	{
		if ( _histogram.size() % 2 == 1 )
			// make number of bins an even number
			_histogram.resize(_histogram.size() + 1, 0);

		// merge bins to reduce size
		for (size_t i=0; i<_histogram.size()/2; i++)
			_histogram[i] = _histogram[2*i] + _histogram[2*i+1];

		_histogram.resize(_histogram.size()/2);

		_interval *= 2;
	}
}

//
// class LinearHistogram
//

template <class Type>
class LinearHistogram : public Histogram<Type> {
public:
    void update(const Type sample);

    LinearHistogram(const string &title, const string &columnLabel, Type interval, size_t maxBins = uMAX, const string &defaultFilename = "") : Histogram<Type>(title, columnLabel, interval, maxBins, defaultFilename) {}
    LinearHistogram(const string &title, const string &columnLabel, Type interval, const string &defaultFilename, size_t maxBins = uMAX)      : Histogram<Type>(title, columnLabel, interval, maxBins, defaultFilename) {}
    LinearHistogram(const string &title, const string &columnLabel, Type interval, const fs::path &defaultFile,   size_t maxBins = uMAX)      : Histogram<Type>(title, columnLabel, interval, maxBins, defaultFile.string()) {}
};

template <class Type>
void LinearHistogram<Type>::update(const Type sample)
{
    Histogram<Type>::update(sample);
	update_histogram(sample);
}

//
// class LogHistogram
//

class LogHistogram : public Histogram<double> {
public:
    enum LogType {B, dB, dBm, dBW};

    void update(const double sample);
    void reset ();

    LogHistogram(const string &title,  const string &columnLabel, double interval, size_t maxBins = uMAX, LogType logType = dB, const string &defaultFilename = "") : Histogram<double>(title, columnLabel, interval, maxBins, defaultFilename),      _logType(logType) {reset();}
    LogHistogram(const string &title,  const string &columnLabel, double interval, const string &defaultFilename, size_t maxBins = uMAX, LogType logType = dB)      : Histogram<double>(title, columnLabel, interval, maxBins, defaultFilename),      _logType(logType) {reset();}
    LogHistogram(const string &title,  const string &columnLabel, double interval, const fs::path &defaultFile,   size_t maxBins = uMAX, LogType logType = dB)      : Histogram<double>(title, columnLabel, interval, maxBins, defaultFile.string()), _logType(logType) {reset();}

private:
    const LogType _logType;

    double _dBSum, _dBSumOfSquares;

    double dBMean()       const {return _count == 0 ? 0.0 : _dBSum / _count;}
    double dBMeanSquare() const {return _count == 0 ? 0.0 : _dBSumOfSquares / _count;}
    double dBStdDev()     const {return sqrt(dBVariance());}
    double dBVariance()   const {return dBMeanSquare() - sqr(dBMean());}
    string dBType()       const;

	void print_footer(ostream &file) const;
};