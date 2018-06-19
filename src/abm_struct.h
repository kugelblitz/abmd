#ifndef DDE_ABM_STRUCT_H
#define DDE_ABM_STRUCT_H

#include "abm.h"

struct _ABM {
  RHS f1;
  RHSD f2;
  int dim;
  double t0;
  double t1;
  double h;
  double *init;
  double *delays;
  int ndelays;
  int abm_order;
  int delays_poly_degree;
  int pointsave_poly_degree;
  double *final_state;
  void *context;
  void (*init_call)(DOUBLE[], void*);
  double *callback_t;
  int (*callback)(double *t, double state[], void *context);
  int **delayed_idxs;
  int *delayed_idxs_lens;
  int *dx_delays_idxs;
  int dx_delays_len;
  char *error;
};

#endif //DDE_ABM_STRUCT_H
