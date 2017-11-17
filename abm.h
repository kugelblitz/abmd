#ifndef DDE_ABM_H
#define DDE_ABM_H

#include "type.h"

DOUBLE *abm(void (*f)(DOUBLE, DOUBLE[], void*, DOUBLE*), DOUBLE, DOUBLE, DOUBLE,
            DOUBLE init[], DOUBLE delays[], int dim, int delays_number, void *context);
void calc_difference(void f(DOUBLE, DOUBLE*, void*, DOUBLE*));

#endif //DDE_ABM_H
