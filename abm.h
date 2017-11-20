#ifndef DDE_ABM_H
#define DDE_ABM_H

#include "type.h"

DOUBLE *abm(void (*f)(DOUBLE t, DOUBLE states[], void *context, DOUBLE *out), DOUBLE, DOUBLE, DOUBLE,
            DOUBLE init[], DOUBLE delays[], int dim, int delays_number, void *context);
#ifdef DEBUG
void calc_difference(void f(DOUBLE, DOUBLE*, void*, DOUBLE*));
#endif

#endif //DDE_ABM_H
