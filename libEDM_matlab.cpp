#include <iostream>

#include <libEDM_library.h>
#include <libEDM_matlab.h>

using std::cerr;
using std::endl;

MATLAB matlab;

cVector convert_mxType (const mxArray *input)
{
    cVector output(mxGetNumberOfElements(input));

    double *inputReal = mxGetPr(input);
    double *inputImag = mxGetPi(input);
    for (size_t index = 0; index < output.size(); index++)
        output[index] = complex<double>(inputReal[index], inputImag[index]);

    return output;
}

mxArray* convert_mxType (const cVector &input)
{
    mxArray *output = mxCreateDoubleMatrix(input.size(), 1, mxCOMPLEX);

    double *outputReal = mxGetPr(output);
    double *outputImag = mxGetPi(output);
    for (size_t index = 0; index < input.size(); index++)
    {
        outputReal[index] = input[index].real();
        outputImag[index] = input[index].imag();
    }

    return output;
}

mxArray* convert_mxType (const cMatrix &input)
{
    mxArray *output = mxCreateDoubleMatrix(input.rows(), input.cols(), mxCOMPLEX);

    double *outputReal = mxGetPr(output);
    double *outputImag = mxGetPi(output);
    for (size_t col = 0; col < input.cols(); col++)
        for (size_t row = 0; row < input.rows(); row++)
        {
            *outputReal++ = input[row][col].real();
            *outputImag++ = input[row][col].imag();
        }

    return output;
}

Engine *MATLAB::engine = NULL;

void MATLAB::openEngine ()
{
    // Open engine if not already open
    if ( !engine )
        engine = engOpen(NULL);

    // Abort if cannot open engine
    if ( !engine )
        error("Cannot open MATLAB engine");

    // Hide MATLAB engine
    engSetVisible(engine, false);
}

MATLAB::~MATLAB ()
{
    if ( engine )
        engClose(engine);
}
