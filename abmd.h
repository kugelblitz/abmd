#ifndef __ABMD_H__
#define __ABMD_H__

#define ABMD_MAX_ORDER 19

#ifdef __MINGW32__
#include <fenv.h>
#define SETENV fesetenv(FE_PC64_ENV)
#else
#define SETENV
#endif

#ifdef USE_LONG_DOUBLE
typedef long double DOUBLE;

#define fmod fmodl
#define sqrt sqrtl
#define pow powl

#else
typedef double DOUBLE;
#endif

typedef struct _ABMD ABMD;

typedef void (*RHS)(DOUBLE x[], double t, DOUBLE *out, void *context);

typedef void (*RHSD)(DOUBLE x[], DOUBLE xs_delayed[], DOUBLE dxs_delayed[],
                     double t, DOUBLE *out, void *context);

ABMD *abmd_create(RHS f, int dim, double t0, double t1, double h, double *init);
void abmd_destroy(ABMD *abm);
int abmd_run(ABMD *);
void abmd_set_order(ABMD *abm, int order);
void abmd_set_delays(ABMD *abm, double *delays, int ndelays);
void abmd_set_delays_poly_degree(ABMD *abm, int deg);
void abmd_set_pointsave_poly_degree(ABMD *abm, int deg);
void abmd_set_f2(ABMD *abm, RHSD f2);
void abmd_set_context(ABMD *abm, void *context);
void abmd_set_init_call(ABMD *abm, void (*init_call)(DOUBLE[], void *));
void abmd_set_callback(ABMD *abm, int (*callback)(double *, double[], void *),
                       double *callback_t);
double *abmd_get_final_state(ABMD *abm);
void abmd_set_delayed_ranges(ABMD *abm, int *ranges, int ranges_len,
                             int delay_idx);
void abmd_set_dx_delays(ABMD *abm, int *idxs, int idxs_len);
const char *abmd_get_last_error(ABMD *abm);

#endif
