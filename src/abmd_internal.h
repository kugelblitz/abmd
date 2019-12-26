#ifndef __ABMD_INTERNAL_H__
#define __ABMD_INTERNAL_H__

#define ABMD_MAX_ORDER 19
#define ABMD_DEFAULT_ORDER 11

#ifdef __MINGW32__
#include <fenv.h>
#define SETENV fesetenv(FE_PC64_ENV)
#else
#define SETENV
#endif

#ifdef ABMD_USE_LONG_DOUBLE
typedef long double ABMD_DOUBLE;
#define fmod fmodl
#define sqrt sqrtl
#define pow powl
#endif

struct _ABMD {
  ABMD_RHS f1;
  ABMD_RHSD f2;
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
  void(*init_call)(ABMD_DOUBLE[], void*);
  double *callback_t;
  int(*callback)(double *t, double state[], void *context);
  int **delayed_idxs;
  int *delayed_idxs_lens;
  int *dx_delays_idxs;
  int dx_delays_len;
  char *error;
};


#endif
