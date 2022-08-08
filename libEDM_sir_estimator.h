#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_matrix.h>

class SIR_Estimator {
public:
	double estimate(const dVector received, const dVector transmitted = dVector()) const;
	double estimate(const dVector received, const bVector transmitted) const;
};