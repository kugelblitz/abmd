#ifndef DDE_RK_H
#define DDE_RK_H

void rk_step(RHS f, double h, double t, DOUBLE *x, int dim, int ndelays,
              void *context, DOUBLE *out, DOUBLE *rhs_out);

#endif
