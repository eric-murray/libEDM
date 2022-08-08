#include <libEDM_sir_estimator.h>

#include <libEDM_library.h>

double SIR_Estimator::estimate(const dVector received, const bVector transmitted) const
{
	dVector _transmitted(transmitted.size());
	for (size_t i=0; i<_transmitted.size(); i++)
		_transmitted[i] = (transmitted[i] ? -1.0 : 1.0);
	return estimate(received, _transmitted);
}

double SIR_Estimator::estimate(const dVector received, const dVector transmitted) const
{
	dVector _transmitted;
	if (transmitted.size() == received.size())
		_transmitted = transmitted;
	else
		_transmitted = received;

	double mean = 0.0;
	for (size_t i=0; i<received.size(); i++)
		mean += received[i]*sgn(_transmitted[i]);
	mean /= received.size();

    double variance = 0.0;
	for (size_t i=0; i<received.size(); i++)
		variance += sqr(received[i]*sgn(_transmitted[i]) - mean);
	variance /= received.size() - 1;

	return sqr(mean) / variance;
}