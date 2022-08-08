#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_library.h>
#include <libEDM_random.h>

Random globalRandom;

radians Random::angle ()
{
    return PI * (2.0*sample() - 1.0);
}

bVector Random::bit_vector (const size_t length)
{
    bVector output;

    for (size_t i=0; i<length; i++)
        output.push_back(bit());

    return output;
}

complex<double> Random::complex_gaussian()
{
    // generates a compex Gaussian deviate with zero mean and unity variance

    // generate random point in unit circle
    double rsq, v1, v2;
    do
    {
        v1 = 2.0*sample() - 1.0;
        v2 = 2.0*sample() - 1.0;
        rsq = sqr(v1) + sqr(v2);
    } while (rsq >= 1.0 || rsq == 0.0);

    // Box-Muller transformation
    double fac = sqrt(-log(rsq)/rsq);

    // return as complex
    return complex<double>(v1 * fac, v2 * fac);
}

double Random::gaussian()
{
    // generates Gaussian deviate with zero mean and unity variance

    if (!iset)
    {
        // compute two independent gaussian samples using complex generator
        complex<double> temp = complex_gaussian();

        // store one sample
        iset = true;
        gset = M_SQRT2 * temp.imag();

        // return the other
        return M_SQRT2 * temp.real();
    }
    else
    {
        // return stored deviate
        iset = false;
        return gset;
    }
}

void Random::lognormal_params( const double mean, const double variance, double &mu, double &sigma )
{
    const double sigma2 = log( 1.0 + variance / mean / mean );
    mu    = std::log( mean ) - 0.5 * sigma2;
    sigma = std::sqrt( sigma2 );
}

size_t Random::poisson (double mean)
{
	// Knuth algorithm
	// Note - computation time is linear with mean. Better algorithms are available

	double p = 0.0;
	size_t k = 0;

	do
	{
		k++;
		p += log(sample());
	}
	while ( p >= -mean );

	return k - 1;
}


void Random::set_seed (const uint32_t seed)
{
    generator.seed(seed);

    uniform01rng = uniform_01<mt19937>(generator);

	iset = false;
}

int Random::uniform (int low, int high)
{
    return static_cast<int>(floor(low + (high - low + 1) * sample()));
}