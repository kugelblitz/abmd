#ifndef DDE_RK_H
#define DDE_RK_H

enum {
  METHOD_DOPRI8 = 1,
  METHOD_RK4 = 2
};

void rk_step(RHS f, double h, double t, DOUBLE *x, int dim, void *context,
             DOUBLE *out, DOUBLE *rhs_out, DOUBLE **memory, int method);

#endif
