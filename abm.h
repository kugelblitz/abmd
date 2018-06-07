#ifndef DDE_ABM_H
#define DDE_ABM_H

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

typedef struct _ABM ABM;

typedef void (*RHS)(DOUBLE x[], double t, DOUBLE *out, void *context);

typedef void (*RHSD)(DOUBLE x[], DOUBLE xs_delayed[], DOUBLE dxs_delayed[],
                     double t, DOUBLE *out, void *context);

int run_abm(ABM *);
ABM *create_abm(RHS f, int dim, double t0, double t1, double h, double *init);
void destroy_abm(ABM *abm);
void set_abm_order(ABM *abm, int order);
void set_delays(ABM *abm, double *delays, int ndelays);
void set_delays_poly_degree(ABM *abm, int deg);
void set_pointsave_poly_degree(ABM *abm, int deg);
void set_f2(ABM *abm, RHSD f2);
void set_context(ABM *abm, void *context);
void set_init_call(ABM *abm, void (*init_call)(DOUBLE[], void*));
void set_callback(ABM *abm, int (*callback)(double *, double[], void*),
                  double *callback_t);
double *get_final_state(ABM *abm);
void set_delayed_ranges(ABM *abm, int *ranges, int ranges_len);
char *get_last_error(ABM *abm);

#endif //DDE_ABM_H
