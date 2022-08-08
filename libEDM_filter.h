#ifndef LIBEDM_FILTER_H
#define LIBEDM_FILTER_H

#include <complex>
#include <deque>

#include <libEDM_matrix.h>

using std::deque;

//
// class Filter
//

template <class State_T, class Coef_T>
class Filter {
public:
    Filter (const Vector<Coef_T> &coefficients) : _coefficients(coefficients), _numCoefficients(coefficients.size()), _state(deque<State_T>(coefficients.size())) {}

    virtual Vector<State_T> filter (const Vector<State_T> &input);

    State_T         operator() (const State_T          input) {return filter(input);}
    Vector<State_T> operator() (const Vector<State_T> &input) {return filter(input);}

    Filter<State_T,Coef_T> operator= (const Filter<State_T,Coef_T> &filter);

    void print(ostream &file = cout, const streamsize precision = 3, const bool final_newline = true) {_coefficients.print(file, precision, final_newline);}

private:
    size_t         _numCoefficients;
    Vector<Coef_T> _coefficients;

    deque<State_T> _state;

    State_T filter (const State_T input);
};

//
// class PulseShapingFilter
//

template <class Type>
class PulseShapingFilter : public Filter<Type,Type> {
public:
    PulseShapingFilter (const Vector<Type> &impulseResponse, const size_t overSamplingFactor) : Filter<Type,Type>(impulseResponse), _overSamplingFactor(overSamplingFactor) {}

    Vector<Type> filter (const Vector<Type> &input) {return Filter<Type,Type>::filter(upsample(input));}

private:
    const size_t _overSamplingFactor;

    Vector<Type> upsample (const Type          input);
    Vector<Type> upsample (const Vector<Type> &input);
};

//
// class RRCFilter
//

class RRCFilter : public PulseShapingFilter<double> {
public:
    RRCFilter (const size_t filterLength, const double rollOffFactor, const size_t overSamplingFactor) : PulseShapingFilter<double>(impulseResponse(filterLength, rollOffFactor, overSamplingFactor), overSamplingFactor) {}

private:
    dVector impulseResponse (const size_t filterLength, const double rollOffFactor, const size_t overSamplingFactor);
};

template <class State_T, class Coef_T>
State_T Filter<State_T,Coef_T>::filter (const State_T input)
{
    // update state with sample
    _state.push_front(input);
    _state.pop_back();

    // compute output
    State_T output = static_cast<State_T>(0.0);
    for (size_t i=0; i<_numCoefficients; i++)
        output += _coefficients[i] * _state[i];

    return output;
}

template <class State_T, class Coef_T>
Vector<State_T> Filter<State_T,Coef_T>::filter (const Vector<State_T> &input)
{
    Vector<State_T> output(input.size());

    for (size_t i=0; i<input.size(); i++)
        output[i] = filter(input[i]);

    return output;
}

template <class State_T, class Coef_T>
Filter<State_T,Coef_T> Filter<State_T,Coef_T>::operator= (const Filter<State_T,Coef_T> &filter)
{
    _numCoefficients = filter._numCoefficients;
    _coefficients    = filter._coefficients;
    _state           = filter._state;

    return (*this);
}

template <class Type>
Vector<Type> PulseShapingFilter<Type>::upsample (const Type input)
{
    Vector<Type> output(_overSamplingFactor);
    output[0] = input;
    return output;
}

template <class Type>
Vector<Type> PulseShapingFilter<Type>::upsample (const Vector<Type> &input)
{
    Vector<Type> output(input.size() * _overSamplingFactor);

    for (size_t i=0; i<input.size(); i++)
        output[i*_overSamplingFactor] = input[i];

    return output;
}

#endif