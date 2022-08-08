#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using std::numeric_limits;
using std::string;
using std::vector;

typedef unsigned long ulong;

// forward declarations
class radians;

template <class Type>
class unit {
public:
	operator double() const {return _value;}

    // explicit method of converting type to double
    double operator() () const {return _value;}

	static Type infinity() {return Type(numeric_limits<double>::infinity());}

	bool operator> (const Type &x) const {return _value > x._value;}
	bool operator< (const Type &x) const {return _value < x._value;}

	Type operator+ (const double x) const {return Type(_value + x);}
	Type operator- (const double x) const {return Type(_value - x);}
	Type operator* (const double x) const {return Type(_value * x);}
	Type operator/ (const double x) const {return Type(_value / x);}

	Type operator+ (const long x) const {return Type(_value + x);}
	Type operator- (const long x) const {return Type(_value - x);}
	Type operator* (const long x) const {return Type(_value * x);}
	Type operator/ (const long x) const {return Type(_value / x);}

	Type operator+ (const size_t x) const {return Type(_value + x);}
	Type operator- (const size_t x) const {return Type(_value - x);}
	Type operator* (const size_t x) const {return Type(_value * x);}
	Type operator/ (const size_t x) const {return Type(_value / x);}

	Type operator- () const {return Type(-_value);}

	Type operator+ (const Type &x) const {return Type(_value + x._value);}
	Type operator- (const Type &x) const {return Type(_value - x._value);}
//	Type operator* (const Type &x) const {return Type(_value * x._value);}
//	Type operator/ (const Type &x) const {return Type(_value / x._value);}

	Type operator+= (const Type &x) {_value += x._value; return Type(_value);}
	Type operator-= (const Type &x) {_value -= x._value; return Type(_value);}
//	Type operator*= (const Type &x) {_value *= x._value; return Type(_value);}
//	Type operator/= (const Type &x) {_value /= x._value; return Type(_value);}

protected:
	double _value;

	explicit unit (const double value) : _value(value) {}
};

class dB : public unit<dB> {
public:
	explicit dB (const double value_dB) : unit(value_dB) {}

	string units() const {return "dB";}
};

class dBi : public dB {
public:
	explicit dBi (const double value_dBi) : dB(value_dBi) {}

	string units() const {return "dBi";}
};

class dBm : public dB {
public:
	static dBm infinity() {return dBm(dB::infinity());}

    explicit dBm (const double value_dBm) : dB(value_dBm) {}

    // explicit declarations of arithmetic operations useful for dBm class
    // those of class unit are not inherited for type dBm
	dBm operator+ (const dB x) const {return dBm(_value + x);}
	dBm operator- (const dB x) const {return dBm(_value - x);}
   	dBm operator- ()           const {return dBm(-_value);}

	string units() const {return "dBm";}
};

class metres : public unit<metres> {
public:
	explicit metres (const double distance_m) : unit(distance_m) {}

	string units() const {return "metres";}
};

class MHz : public unit<MHz> {
public:
	explicit MHz (const double frequency_MHz) : unit(frequency_MHz) {}

	string units() const {return "MHz";}
};

class degrees : public unit<degrees> {
    friend class unit<degrees>;

public:
    degrees ()                     : unit(0.0)          {}
    degrees (const degrees &angle) : unit(bound(angle)) {}

    degrees (const radians &angle);

    string units() const {return "degrees";}

private:
    degrees (const double angle) : unit(bound(angle)) {}

    double bound (double value) const;
};

class radians : public unit<radians> {
    friend class unit<radians>;
    template <class Type> friend Type atan (const double x);
    template <class Type> friend Type atan2(const double x, const double y);
    template <class Type> friend Type acos (const double x);
    template <class Type> friend Type asin (const double x);

public:
    radians ()                     : unit(0.0)          {}
    radians (const radians &angle) : unit(bound(angle)) {}

    radians (const degrees &angle);

    string units() const {return "radians";}

private:
    radians (const double angle) : unit(bound(angle)) {}

    double bound (double value) const;
};

//
// PointerVector
//

template <class Type>
class PointerVector : protected vector<Type*> {
public:
    using vector<Type*>::at;
	using vector<Type*>::back;
	using vector<Type*>::begin;
    using vector<Type*>::clear;
    using vector<Type*>::end;
    using vector<Type*>::front;
    using vector<Type*>::operator[];
    using vector<Type*>::pop_back;
    using vector<Type*>::push_back;
    using vector<Type*>::size;

    PointerVector(const bool deleteObjectsWhenFinished = true) : _deleteObjectsWhenFinished(deleteObjectsWhenFinished) {}
    ~PointerVector() {if (_deleteObjectsWhenFinished) clearAndDelete();}

    void clearAndDelete();

private:
    const bool _deleteObjectsWhenFinished;
};

template <class Type>
void PointerVector<Type>::clearAndDelete()
{
    while (!empty())
    {
        Type *object = back();
        pop_back();
        delete object;
    }
}

const degrees oneDegree = degrees() + 1.0;
const radians PI        = radians() + M_PI;