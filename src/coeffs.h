#ifndef ADAMS_COEFFS_H
#define ADAMS_COEFFS_H

#include "abm.h"

void get_corrector_coeffs(int, DOUBLE *);
void get_predictor_coeffs(int, DOUBLE *);

extern DOUBLE PREDICTOR_COEFFS[19];

#endif //ADAMS_COEFFS_H
