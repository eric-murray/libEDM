#pragma once

#include <complex>
#include <cstdlib>
#include <ctime>
#include <vector>

#include <boost/cstdint.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>

#include <libEDM_types.h>
#include <libEDM_matrix.h>

using std::complex;
using std::time;

using boost::uint32_t;
using boost::mt19937;
using boost::uniform_01;
using boost::uniform_int;

class Random {
public:
    // constructor
    Random (uint32_t seed = time(0)) : generator(seed), iset(false), uniform01rng(generator) {}

    // generates random number between (0, 1]
    double sample () { return uniform01rng(); }

    // generates a random angle between (-PI, PI]
    radians angle ();

    // generates uniform random number between (low, high)
    int uniform (int low, int high);

    // generates a Gaussian random number with specified mean and variance
    double gaussian (double mean, double stddev) { return mean + (stddev == 0.0 ? 0.0 : stddev * gaussian()); }
    double gaussian ();

    // generates a complex Gaussian random number with specified mean and variance
    complex<double> complex_gaussian (const complex<double> &mean, double stddev) { return mean + (stddev == 0.0 ? 0.0 : stddev * complex_gaussian()); }
    complex<double> complex_gaussian ();

    // lognormal distribution
    double lognormal( const double mu, const double sigma ) {return std::exp( gaussian(mu, sigma) );}

    void lognormal_params( const double mean, const double variance, double &mu, double &sigma );

	// generates a Poisson random number
	size_t poisson (double mean);

    // generates a Rayleigh random number
    double rayleigh (double stddev = 1.0) {return norm(complex_gaussian(0.0, stddev));}

    // generates an Exponential random number
    double exponential (const double mu) {return -mu * std::log( 1.0 - sample() );}

    // Pareto distribution
    double pareto ( const double alpha, const double min ) {return min * std::pow( sample(), -1.0 / alpha );}

    double pareto_min ( const double k, const double mean ) {return ( k - 1.0 ) * mean / k;}

    // generates a random bit
    bool bit () { return sample() < 0.5; }

    // generates a vector of random bits
    bVector bit_vector (const size_t length);

    // set seed of random generator
    void set_seed (const uint32_t seed);

private:
    mt19937             generator;
    uniform_01<mt19937> uniform01rng;

    bool   iset;
    double gset;
};

extern Random globalRandom;