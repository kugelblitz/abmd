#ifndef DDE_RK_H
#define DDE_RK_H

void rk_step(RHS f, double h, double t, DOUBLE *x, int dim, int ndelays,
             int *delayed_idxs, int delayed_idxs_len, void *context,
             DOUBLE *out, DOUBLE *rhs_out, DOUBLE **memory);

#endif
