#ifndef DDE_RK_H
#define DDE_RK_H

enum {
  METHOD_DOPRI8 = 1,
  METHOD_RK4 = 2
};

void rk_step(ABMD_RHS f, double h, double t, ABMD_DOUBLE *x, int dim, void *context,
             ABMD_DOUBLE *out, ABMD_DOUBLE *rhs_out, ABMD_DOUBLE **memory, int method);

#endif
