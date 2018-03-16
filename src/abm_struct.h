#ifndef DDE_ABM_STRUCT_H
#define DDE_ABM_STRUCT_H

#include "abm.h"

struct _ABM {
  void (*f1)(DOUBLE *, double, DOUBLE *, void *);
  void (*f2)(DOUBLE *, double, DOUBLE *, void *);
  int dim;
  double t0;
  double t1;
  double h;
  double *init;
  double *delays;
  int ndelays;
  int _user_set_delays;
  int abm_order;
  int interpolation_order;
  int extrapolation_order;
  double *final_state;
  void *context;
  void (*init_call)(DOUBLE[], void*);
  double *callback_t;
  int (*callback)(double *t, double state[], void *context);
};

#endif //DDE_ABM_STRUCT_H
