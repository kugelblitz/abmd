#ifndef DDE_ABM_STRUCT_H
#define DDE_ABM_STRUCT_H

#include "abm.h"

struct _ABM {
  RHS1 f1;
  RHS2 f2;
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
  int *delayed_idxs;
  int delayed_idxs_len;
};

#endif //DDE_ABM_STRUCT_H
