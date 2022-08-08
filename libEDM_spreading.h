#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_matrix.h>

template <class VectorType>
class Spreader {
public:
	Spreader (const dMatrix &spreadingCodes) : spreadingCodes(spreadingCodes) {}

	void       spread   (const VectorType &input, const size_t spreadingCodeId, VectorType &output) const;
    VectorType spread   (const VectorType &input, const size_t spreadingCodeId) const;
    void       despread (const VectorType &input, const size_t spreadingCodeId, VectorType &output) const;
    VectorType despread (const VectorType &input, const size_t spreadingCodeId) const;

    size_t SF() const {return spreadingCodes.size();}

private:
	const dMatrix spreadingCodes;
};

template <class VectorType>
void Spreader<VectorType>::spread (const VectorType &input, const size_t spreadingCodeId, VectorType &output) const
{
    output.resize(0);
    for (size_t i=0; i<input.size(); i++)
        output.ins(output.size(), spreadingCodes[spreadingCodeId] * input[i]);
}

template <class VectorType>
VectorType Spreader<VectorType>::spread (const VectorType &input, const size_t spreadingCodeId) const
{
    VectorType output;
    spread(input, spreadingCodeId, output);
    return output;
}

template <class VectorType>
void Spreader<VectorType>::despread (const VectorType &input, const size_t spreadingCodeId, VectorType &output) const
{
    output.clear();
    output.resize(input.size() / SF());

    for (size_t i=0; i<output.size(); i++)
    {
        for (size_t j=0; j<SF(); j++)
            output[i] += input[i*SF() + j] * spreadingCodes[spreadingCodeId][j];
        output[i] /= SF();
    }
}

template <class VectorType>
VectorType Spreader<VectorType>::despread (const VectorType &input, const size_t spreadingCodeId) const
{
    VectorType output;
    despread(input, spreadingCodeId, output);
    return output;
}

dMatrix wcdma_spreading_codes(const size_t SF);