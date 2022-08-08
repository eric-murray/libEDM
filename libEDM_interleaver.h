#pragma once

#include <assert.h>

#include <libEDM_matrix.h>

template <class Type>
class SequenceInterleaver {
public:
    Vector<Type> interleave (const Vector<Type> &input) const;
    void         interleave (const Vector<Type> &input, Vector<Type> &output) const;

    Vector<Type> deinterleave (const Vector<Type> &input) const;
    void         deinterleave (const Vector<Type> &input, Vector<Type> &output) const;

    uVector get_interleaver_sequence(void) const {return interleaver_sequence;}
    size_t  interleaver_depth(void)        const {return interleaver_sequence.size();}

    void set_interleaver_sequence(uVector interleaver_sequence) {this->interleaver_sequence = interleaver_sequence;}

    SequenceInterleaver(const uVector interleaver_sequence) : interleaver_sequence(interleaver_sequence) {}

private:
    uVector interleaver_sequence;
};

typedef SequenceInterleaver<bool>   bSequenceInterleaver;
typedef SequenceInterleaver<double> dSequenceInterleaver;

template <class Type>
void SequenceInterleaver<Type>::interleave(const Vector<Type> &input, Vector<Type> &output) const
{
    size_t steps         = static_cast<size_t>(ceil(static_cast<double>(input.size()) / interleaver_depth()));
    size_t output_length = steps * interleaver_depth();

    assert( output_length == input.size() );

    output.resize(output_length);
    for (size_t s=0; s<steps; s++)
        for (size_t i=0; i<interleaver_depth(); i++)
            output[s*interleaver_depth() + i] = input[s*interleaver_depth() + interleaver_sequence[i]];
}

template <class Type>
Vector<Type> SequenceInterleaver<Type>::interleave(const Vector<Type> &input) const
{
    Vector<Type> output;
    interleave(input, output);
    return output;
}

template <class Type>
void SequenceInterleaver<Type>::deinterleave(const Vector<Type> &input, Vector<Type> &output) const
{
    size_t steps         = static_cast<size_t>(ceil(divide(input.size(), interleaver_depth())));
    size_t output_length = steps * interleaver_depth();

    assert ( output_length == input.size() );

    output.resize(output_length);
	for (size_t s=0; s<steps; s++)
	    for (size_t i=0; i<interleaver_depth(); i++)
	        output[s*interleaver_depth() + interleaver_sequence[i]] = input[s*interleaver_depth() + i];
}

template <class Type>
Vector<Type> SequenceInterleaver<Type>::deinterleave(const Vector<Type> &input) const
{
    Vector<Type> output;
    deinterleave(input, output);
    return output;
}