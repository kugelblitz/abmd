#ifndef __ABMD_H__
#define __ABMD_H__

typedef struct _ABMD ABMD;

#ifdef ABMD_USE_LONG_DOUBLE
typedef long double ABMD_DOULBgE;
#else
typedef double ABMD_DOUBLE;
#endif


typedef void (*ABMD_RHS)(ABMD_DOUBLE x[], double t, ABMD_DOUBLE *out, void *context);

typedef void (*ABMD_RHSD)(ABMD_DOUBLE x[], ABMD_DOUBLE xs_delayed[], ABMD_DOUBLE dxs_delayed[],
                     double t, ABMD_DOUBLE *out, void *context);

ABMD *abmd_create(ABMD_RHS f, int dim, double t0, double t1, double h, double *init);
void abmd_destroy(ABMD *abm);
int abmd_run(ABMD *);
void abmd_set_order(ABMD *abm, int order);
void abmd_set_delays(ABMD *abm, double *delays, int ndelays);
void abmd_set_delays_poly_degree(ABMD *abm, int deg);
void abmd_set_pointsave_poly_degree(ABMD *abm, int deg);
void abmd_set_f2(ABMD *abm, ABMD_RHSD f2);
void abmd_set_context(ABMD *abm, void *context);
void abmd_set_init_call(ABMD *abm, void (*init_call)(ABMD_DOUBLE[], void *));
void abmd_set_callback(ABMD *abm, int (*callback)(double *, double[], void *),
                       double *callback_t);
double *abmd_get_final_state(ABMD *abm);
void abmd_set_delayed_ranges(ABMD *abm, int *ranges, int ranges_len,
                             int delay_idx);
void abmd_set_dx_delays(ABMD *abm, int *idxs, int idxs_len);
const char *abmd_get_last_error(ABMD *abm);

#endif
