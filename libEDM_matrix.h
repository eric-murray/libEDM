#pragma once

#include <cassert>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

#include <libEDM_types.h>

using std::abs;
using std::complex;
using std::conj;
using std::cout;
using std::endl;
using std::fixed;
using std::nouppercase;
using std::ostream;
using std::scientific;
using std::setprecision;
using std::setw;
using std::streamsize;
using std::uppercase;
using std::vector;

template <class Type>
class Vector : public vector<Type> {
public:
	//using vector<Type>::assign;
	//using vector<Type>::begin;
	//using vector<Type>::const_iterator;
	//using vector<Type>::clear;
	//using vector<Type>::end;
	//using vector<Type>::operator =;
	//using vector<Type>::operator [];
	//using vector<Type>::push_back;
	//using vector<Type>::reserve;
	//using vector<Type>::resize;

    Vector(const size_t num = 0, const Type &val = static_cast<Type>(0)) : vector<Type>(num, val) {}
    Vector(const Vector<Type> &vec)                                      : vector<Type>(vec)      {}
    Vector(const Type *array, const size_t num);

    void print (ostream &file = cout, const streamsize precision = 3, const bool final_newline = true) const;

    Type   sum           () const;
    Type   sum_of_modulus() const;
    Type   sum_of_squares() const;
	double mean          () const {return static_cast<double>(sum())            / size();}
	double mean_square   () const {return static_cast<double>(sum_of_squares()) / size();}
    double variance      () const;

    void zeros();
    void ones ();
    void index();
    void ramp (const Type start, const Type end);
    void ramp (const Type start, const Type step, const Type range, const Type max = numeric_limits<Type>::infinity());

    void add      (const Vector<Type> &x);
    void subtract (const Vector<Type> &x);
    void multiply (const Vector<Type> &x);
    void or       (const Vector<Type> &x);
    void xor      (const Vector<Type> &x);
    void and      (const Vector<Type> &x);
    void not      ();
    void conj     ();

    void multiply (const Type &x);

    Vector<Type> operator+ (const Vector<Type> &x) const;
    Vector<Type> operator+=(const Vector<Type> &x);
    Vector<Type> operator- (const Vector<Type> &x) const;
    Vector<Type> operator-=(const Vector<Type> &x);
    Vector<Type> operator* (const Vector<Type> &x) const;
    Vector<Type> operator*=(const Vector<Type> &x);
    Vector<Type> operator| (const Vector<Type> &x) const;
    Vector<Type> operator|=(const Vector<Type> &x);
    Vector<Type> operator^ (const Vector<Type> &x) const;
    Vector<Type> operator^=(const Vector<Type> &x);
    Vector<Type> operator& (const Vector<Type> &x) const;
    Vector<Type> operator&=(const Vector<Type> &x);
    Vector<Type> operator! ()                      const;

	Vector<Type>           operator/ (const Type &x) const;
    Vector<Type>           operator/=(const Type &x);
    Vector<Type>           operator* (const Type &x) const;
    Vector<Type>           operator*=(const Type &x);
    Vector<complex<Type> > operator* (const complex<Type> &x) const;

    Vector<Type> mid  (size_t start, size_t nr) const;
    Vector<Type> right(size_t nr)               const;
    Vector<Type> left (size_t nr)               const;

    void replace_mid(const size_t pos, const Vector<Type> &x);
    void del        (const size_t index);
    void ins        (const size_t index, const Type &in);
    void ins        (const size_t index, const Vector<Type> &in);

};

typedef Vector<double>           dVector;
typedef Vector<int>              iVector;
typedef Vector<bool>             bVector;
typedef Vector<complex<double> > cVector;
typedef Vector<size_t>           uVector;

template <class Type>
class Matrix : public vector<Vector<Type> > {
public:
    size_t rows() const {return size();}
    size_t cols() const {return at(0).size();}

    Matrix<Type> (size_t num = 0);
    Matrix<Type> (size_t num, Type val)                             {set_size(num, num, val);}
    Matrix<Type> (size_t rows, size_t columns, Type val)            {set_size(rows, columns, val);}
    Matrix<Type> (const Matrix<Type> &x) : vector<Vector<Type> >(x) {};

    void add      (const Matrix<Type> &x);
    void subtract (const Matrix<Type> &x);
    void invert   ();
    void multiply (const Matrix<Type> &x);

    Matrix<Type> operator+  (const Matrix<Type> &x) const;
    Matrix<Type> operator+= (const Matrix<Type> &x);
    Matrix<Type> operator-  (const Matrix<Type> &x) const;
    Matrix<Type> operator-= (const Matrix<Type> &x);
    Matrix<Type> operator*  (const Matrix<Type> &x) const;
    Matrix<Type> operator*  (const Type &x)         const;
    Matrix<Type> operator*= (const Matrix<Type> &x);
    Matrix<Type> operator!  ()                      const;
    Matrix<Type> inverse    ()                      const;
	Matrix<Type> transpose  ()                      const;

    Vector<Type> operator*  (const Vector<Type> &x) const;

    void zeros    ();
    void set_size (const size_t rows, const size_t columns, const Type val = static_cast<Type>(0));
    void print    (ostream &file = cout, const streamsize precision = 3) const;

};

typedef Matrix<bool>             bMatrix;
typedef Matrix<double>           dMatrix;
typedef Matrix<int>              iMatrix;
typedef Matrix<size_t>           uMatrix;
typedef Matrix<complex<double> > cMatrix;

template <class Type>
class Cubrix : public vector<Matrix<Type> > {
public:
    Cubrix<Type> (const size_t dim1, const size_t dim2, const size_t dim3, const Type val = static_cast<Type>(0)) {set_size(dim1, dim2, dim3, val);}

private:
    void set_size (const size_t dim1, const size_t dim2, const size_t dim3, const Type val = static_cast<Type>(0));
};

typedef Cubrix<bool>   bCubrix;
typedef Cubrix<double> dCubrix;

template<class Type>
Type inner_product(const Vector<Type> &x, const Vector<Type> &y)
{
	assert( x.size() == y.size() );

	// Compute vector dot product
	Type output = 0;
	for (size_t i=0; i<x.size(); i++)
		output += x[i] * y[i];

	return output;
}

template<class Type>
Matrix<Type> outer_product(const Vector<Type> &x, const Vector<Type> &y)
{
	assert( x.size() == y.size() );

	// Compute outer product
	Matrix<Type> output(x.size(), 0);
	for (size_t i=0; i<x.size(); i++)
		for (size_t j=0; j<x.size(); j++)
			output[i][j] = x[i] * y[j];

	return output;
}


template <class Type>
Vector<Type>::Vector(const Type *array, const size_t num) : vector<Type>(num)
{
    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = array[i];
}

template<class Type>
void Vector<Type>::print (ostream &file, const streamsize precision, const bool final_newline) const
{
    streamsize width;
    if (precision == 0)
    {
        file << fixed;
        width = 0;
    }
    else
    {
        file << scientific << uppercase;
        width = 10;
    }

    file << setprecision(precision);
    for (size_t i=0; i<this->size(); i++)
        file << setw(width) << this->at(i) << " ";
    file << fixed << nouppercase;

    if (final_newline)
        file << endl;
}

template <class Type>
Type Vector<Type>::sum() const
{
    Type output = static_cast<Type>(0);
    for (size_t i=0; i<this->size(); i++)
        output += this->at(i);
    return output;
}

template <class Type>
Type Vector<Type>::sum_of_modulus() const
{
    Type output = static_cast<Type>(0);
    for (size_t i=0; i<this->size(); i++)
		if (this->at(i) > static_cast<Type>(0))
            output += this->at(i);
		else
			output -= this->at(i);
    return output;
}

template <class Type>
Type Vector<Type>::sum_of_squares() const
{
    Type output = static_cast<Type>(0);
    for (size_t i=0; i<this->size(); i++)
        output += sqr(this->at(i));
    return output;
}

template <class Type>
double Vector<Type>::variance() const
{
    double output = 0.0;
    for (size_t i=0; i<size(); i++)
        output += sqr(this->at(i) - mean());
    output /= (size() - 1);
    return output;
}

template<class Type>
void Vector<Type>::zeros()
{
    for (size_t i=0; i<this->size(); i++)
        this->at(i) = static_cast<Type>(0);
}

template<class Type>
void Vector<Type>::ones()
{
    for (size_t i=0; i<this->size(); i++)
        this->at(i) = static_cast<Type>(1);
}

template<class Type>
void Vector<Type>::index()
{
    for (size_t i=0; i<this->size(); i++)
        this->at(i) = static_cast<Type>(i);
}

template <class Type>
void Vector<Type>::ramp (const Type start, const Type end)
{
    Type step = (end - start) / static_cast<Type>(size()-1);
    ramp(start, step, 0.0);
}

template <class Type>
void Vector<Type>::ramp (const Type start, const Type step, const Type range, const Type max)
{
    for (size_t i=0; i<size(); i++)
    {
        Type value = start + i*step;
        while (value > max)
            value -= range;
        this->at(i) = value;
    } 
}

template <class Type>
void Vector<Type>::add (const Vector<Type> &x)
{
    // check that vectors are the same size
    assert ( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = (*this)[i] + x[i];
}

template <class Type>
Vector<Type> Vector<Type>::operator+ (const Vector<Type> &x) const
{
    Vector<Type> output = *this;
    output.add(x);
    return output;
}

template <class Type>
Vector<Type> Vector<Type>::operator+= (const Vector<Type> &x)
{
    this->add(x);
    return *this;
}

template <class Type>
void Vector<Type>::subtract (const Vector<Type> &x)
{
    // check that vectors are the same size
    assert ( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] -= x[i];
}

template <class Type>
Vector<Type> Vector<Type>::operator- (const Vector<Type> &x) const
{
    Vector<Type> output = *this;
    output.subtract(x);
    return output;
}

template <class Type>
Vector<Type> Vector<Type>::operator-= (const Vector<Type> &x)
{
    this->subtract(x);
    return *this;
}

template <class Type>
void Vector<Type>::multiply (const Vector<Type> &x)
{
    // check that vectors are the same size
    assert( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = (*this)[i] * x[i];
}

template <class Type>
Vector<Type> Vector<Type>::operator* (const Vector<Type> &x) const
{
    Vector<Type> output = *this;
    output.multiply(x);
    return output;
}

template <class Type>
Vector<Type> Vector<Type>::operator*= (const Vector<Type> &x)
{
    this->multiply(x);
    return *this;
}

template <class Type>
void Vector<Type>::multiply (const Type &x)
{
    for (size_t i=0; i<this->size(); i++)
        (*this)[i] *= x;
}

template <class Type>
Vector<Type> Vector<Type>::operator* (const Type &x) const
{
    Vector<Type> output = *this;
    output.multiply(x);
    return output;
}

template <class Type>
Vector<Type> Vector<Type>::operator*= (const Type &x)
{
    this->multiply(x);
    return *this;
}

template <class Type>
Vector<Type> Vector<Type>::operator/ (const Type &x) const
{
    return operator*(static_cast<Type>(1) / x);
}

template <class Type>
Vector<Type> Vector<Type>::operator/= (const Type &x)
{
    return operator*=(static_cast<Type>(1) / x);
}

template <class Type>
Vector<complex<Type> > Vector<Type>::operator* (const complex<Type> &x) const
{
    Vector<complex<Type> > output;
    for (size_t i=0; i<this->size(); i++)
        output.push_back(x * this->at(i));

    return output;
}

template <class Type>
void Vector<Type>::xor (const Vector<Type> &x)
{
    // check that vectors are the same size
    assert( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = (*this)[i] ^ x[i];
}

template <class Type>
Vector<Type> Vector<Type>::operator^ (const Vector<Type> &x) const
{
    Vector<Type> output = *this;
    output.xor(x);
    return output;
}

template <class Type>
Vector<Type> Vector<Type>::operator^= (const Vector<Type> &x)
{
    this->xor(x);
    return *this;
}

template <class Type>
void Vector<Type>::or (const Vector<Type> &x)
{
    // check that vectors are the same size
    assert( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = (*this)[i] | x[i];
}

template <class Type>
Vector<Type> Vector<Type>::operator| (const Vector<Type> &x) const
{
    Vector<Type> output = *this;
    output.or(x);
    return output;
}

template <class Type>
Vector<Type> Vector<Type>::operator|= (const Vector<Type> &x)
{
    this->or(x);
    return *this;
}

template <class Type>
void Vector<Type>::and (const Vector<Type> &x)
{
    // check that vectors are the same size
    assert( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = (*this)[i] & x[i];
}

template <class Type>
Vector<Type> Vector<Type>::operator& (const Vector<Type> &x) const
{
    Vector<Type> output = *this;
    output.and(x);
    return output;
}

template <class Type>
Vector<Type> Vector<Type>::operator&= (const Vector<Type> &x)
{
    this->and(x);
    return *this;
}

template <class Type>
void Vector<Type>::not ()
{
    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = !(*this)[i];
}

template <class Type>
Vector<Type> Vector<Type>::operator! () const
{
    Vector<Type> output = *this;
    output.not();
    return output;
}

template <class Type>
void Vector<Type>::conj()
{
    for (size_t i=0; i<this->size(); i++)
        (*this)[i] = std::conj((*this)[i]);
}

template <class Type>
Vector<Type> Vector<Type>::mid(size_t start, size_t nr) const
{
    // check that there are nr elements to return
    assert( (start + nr) <= this->size() );

    Vector<Type> output;
    for (size_t i=start; i<(start+nr); i++)
        output.push_back(this->at(i));

    return output;
}

template <class Type>
Vector<Type> Vector<Type>::right(size_t nr) const
{
    return mid(this->size() - nr, nr);
}

template <class Type>
Vector<Type> Vector<Type>::left(size_t nr) const
{
    return mid(0, nr);
}

template <class Type>
void Vector<Type>::replace_mid(const size_t pos, const Vector<Type> &x)
{
    // check that elements to be replaced exist
    assert( (pos >= 0) && (pos + x.size() <= this->size()) );

    for (size_t i=0; i<x.size(); i++)
        this->at(pos+i) = x[i];
}

template <class Type>
void Vector<Type>::del(const size_t index)
{
    // check that element to be deleted exists
    assert( (index >= 0) && (index < this->size()) );

    this->erase(this->begin() + index);
}

template <class Type>
void Vector<Type>::ins(const size_t index, const Type &in)
{
    // check that index is valid
    assert( (index >= 0) && (index <= this->size()) );

    this->insert(this->begin() + index, in);
}

template <class Type>
void Vector<Type>::ins(const size_t index, const Vector<Type> &in)
{
    // check that index is valid
    assert( (index >= 0) && (index <= this->size()) );

    this->insert(this->begin() + index, in.begin(), in.end());
}

template <class Type>
Matrix<Type>::Matrix (size_t num)
{
    // if no val specified, initialise to identity Matrix<Type>
    set_size(num, num);

    for (size_t i=0; i<num; i++)
        (*this)[i][i] = static_cast<Type>(1.0);
}

template <class Type>
void Matrix<Type>::print(ostream &file, const streamsize precision) const
{
    for (size_t i=0; i<this->size(); i++)
        this->at(i).print(file, precision, true);
    file << endl;
}

template <class Type>
void Matrix<Type>::zeros()
{
    for (size_t i=0; i<this->size(); i++)
        (*this)[i].zeros();
}

template <class Type>
void Matrix<Type>::set_size (const size_t rows, const size_t columns, const Type val)
{
    this->clear();
    this->resize(rows);
    for (size_t i=0; i<this->size(); i++)
        (*this)[i].resize(columns, val);
}

template <class Type>
void Matrix<Type>::add (const Matrix<Type> &x)
{
    // check that vectors are the same size
    assert( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] += x[i];
}

template <class Type>
Matrix<Type> Matrix<Type>::operator+ (const Matrix<Type> &x) const
{
    Matrix<Type> output = *this;
    output.add(x);
    return output;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator+= (const Matrix<Type> &x)
{
    this->add(x);
    return *this;
}

template <class Type>
void Matrix<Type>::subtract (const Matrix<Type> &x)
{
    // check that vectors are the same size
    assert( this->size() == x.size() );

    for (size_t i=0; i<this->size(); i++)
        (*this)[i] -= x[i];
}

template <class Type>
Matrix<Type> Matrix<Type>::operator- (const Matrix<Type> &x) const
{
    Matrix<Type> output = *this;
    output.subtract(x);
    return output;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator-= (const Matrix<Type> &x)
{
    this->subtract(x);
    return *this;
}

template <class Type>
Matrix<Type> Matrix<Type>::inverse() const
{
    Matrix<Type> output = *this;
    output.invert();
    return output;
}

template <class Type>
Matrix<Type> Matrix<Type>::transpose() const
{
	Matrix<Type> output(this->cols(), this->rows(), 0.0);

	for (size_t row=0; row<output.rows(); row++)
		for (size_t col=0; col<output.cols(); col++)
			output[row][col] = (*this)[col][row];

	return output;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator! () const
{
    return inverse();
}

template <class Type>
void Matrix<Type>::invert()
{
    // Matrix<Type> inversion by Gaussian elimination without pivoting

    // store Matrix<Type> size in local variable
    const size_t n = this->size();

    for (size_t i=0; i<n; i++)
    {
        // assert that pivot value is not zero
        assert( this->at(i)[i] != 0.0 );

        // invert pivot value
        (*this)[i][i] = 1.0 / (*this)[i][i];

        // multiply other elements in same column by inverted pivot value
        for (size_t j=0; j<n; j++)
            if (i != j)
                (*this)[i][j] *= (*this)[i][i];

        // iterate over other columns
        for (size_t k=0; k<n; k++)
            if (k != i)
            {
                // iterate over rows
                for (size_t j=0; j<n; j++)
                    if (j != i)
                        // element is not in same row as pivot value
                        (*this)[k][j] -= (*this)[k][i] * (*this)[i][j];

                // negate element in same row as pivot, and multiply by inverted pivot value
                // note that this operation must be performed after output[k][i] has been used to modify other elements
                (*this)[k][i] = -(*this)[k][i] * (*this)[i][i];
            }
    }
}


template <class Type>
void Matrix<Type>::multiply (const Matrix<Type> &x)
{
    // multiplies this by x, modifying this in the process

    // check that matrices are of equal dimension
    assert( this->size() == x.size() );

    // create temporary Matrix<Type> to store contents of multiplication
    Matrix<Type> temp(this->size(), static_cast<Type>(0));

    // multiply temp by x and copy to this
    for (size_t i=0; i<this->size(); i++)
        for (size_t j=0; j<this->size(); j++)
            for (size_t k=0; k<this->size(); k++)
                temp[i][j] += (*this)[i][k] * x[k][j];

    // copy result to this
    for (size_t i=0; i<this->size(); i++)
        for (size_t j=0; j<this->size(); j++)
            (*this)[i][j] = temp[i][j];
}

template <class Type>
Vector<Type> Matrix<Type>::operator* (const Vector<Type> &x) const
{
    // multiplies Matrix<Type> by column vector x and returns result

    // check that inner dimensions agree
    assert( this->size() == x.size() );

    // define output vector
    Vector<Type> output(x.size(), static_cast<Type>(0.0));

    for (size_t i=0; i<output.size(); i++)
        for (size_t k=0; k<x.size(); k++)
            output[i] += (*this)[i][k] * x[k];

    return output;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator* (const Type &x) const
{
    Matrix<Type> output = *this;
	for (size_t i=0; i<this->size(); i++)
		output[i] *= x;
    return output;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator* (const Matrix<Type> &x) const
{
    Matrix<Type> output = *this;
    output.multiply(x);
    return output;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator*= (const Matrix<Type> &x)
{
    this->multiply(x);
    return *this;
}

template <class Type>
void Cubrix<Type>::set_size (const size_t dim1, const size_t dim2, const size_t dim3, const Type val = static_cast<Type>(0))
{
    clear();
    resize(dim1);
    for (size_t i=0; i<this->size(); i++)
    {
        (*this)[i].resize(dim2);
         for (size_t j=0; j<(*this)[i].size(); j++)
             (*this)[i][j].resize(dim3, val);
    }
}