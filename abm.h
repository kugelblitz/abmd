#ifndef DDE_ABM_H
#define DDE_ABM_H

#ifdef USE_LONG_DOUBLE
typedef long double DOUBLE;

#define fmod fmodl
#define sqrt sqrtl
#define pow powl
#define ceil ceill

#else
typedef double DOUBLE;
#endif

typedef struct _ABM ABM;

void run_abm(ABM *);
ABM *create_abm(void (*f)(double, DOUBLE[], void*, DOUBLE*), int dim,
                double t0, double t1, double h, double *init);
void destroy_abm(ABM *abm);
void set_abm_order(ABM *abm, int order);
void set_delays(ABM *abm, double *delays, int ndelays);
void set_interpolation_order(ABM *abm, int order);
void set_extrapolation_order(ABM *abm, int order);
void set_f2(ABM *abm, void (*f2)(double, DOUBLE[], void*, DOUBLE*));
void set_context(ABM *abm, void *context);
void set_init_call(ABM *abm, void (*init_call)(DOUBLE[], void*));
void set_callback(ABM *abm, int (*callback)(double *, double[], void*),
                  double *callback_t);
double *get_final_state(ABM *abm);

#endif //DDE_ABM_H
