#pragma once

#include <engine.h>


#include <libEDM_matrix.h>

class MATLAB {
public:
    static Engine *engine;

    ~MATLAB();

    static void openEngine();
};

extern MATLAB matlab;

cVector  convert_mxType (const mxArray *input);
mxArray* convert_mxType (const cVector &input);
mxArray* convert_mxType (const cMatrix &input);