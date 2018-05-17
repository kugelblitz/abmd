#include <stdlib.h>
#include <string.h>

#include "abm.h"
#include "abm_struct.h"

#define ABM_ORDER 11

ABM *create_abm(RHS1 f, int dim, double t0, double t1, double h, double *init) {
  ABM *abm = (ABM *) malloc(sizeof(ABM));
  double *delays = (double *) malloc(sizeof(double));
  delays[0] = 0;
  double *final_state = (double *) malloc(sizeof(double) * dim);
  *abm = (ABM) {
          .f1=f,
          .f2=NULL,
          .dim=dim,
          .t0=t0,
          .t1=t1,
          .h=h,
          .init=init,
          .delays=delays,
          .ndelays=1,
          .abm_order=ABM_ORDER,
          .interpolation_order=4,
          .extrapolation_order=1,
          .final_state=final_state,
          .context=NULL,
          .init_call=NULL,
          .callback_t=NULL,
          .callback=NULL,
          .delayed_idxs=NULL,
          .delayed_idxs_len = dim
  };
  return abm;
}

void destroy_abm(ABM *abm) {
  free(abm->delays);
  free(abm->final_state);
  free(abm->delayed_idxs);
  free(abm);
}

void set_abm_order(ABM *abm, int order) {
  abm->abm_order = order;
}

void set_delays(ABM *abm, double *delays, int ndelays) {
  if (abm->delays != NULL) {
    free(abm->delays);
    abm->delays = NULL;
  }
  abm->ndelays = ndelays;
  if (ndelays > 0) {
    abm->delays = malloc(sizeof(double) * ndelays);
    memcpy(abm->delays, delays, sizeof(double) * ndelays);
  }
}

void set_interpolation_order(ABM *abm, int order) {
  abm->interpolation_order = order;
}

void set_extrapolation_order(ABM *abm, int order) {
  abm->extrapolation_order = order;
}

void set_f2(ABM *abm, RHS2 f2) {
  abm->f2 = f2;
}

void set_context(ABM *abm, void *context) {
  abm->context = context;
}

void set_init_call(ABM *abm, void (*init_call)(DOUBLE[], void*)) {
  abm->init_call = init_call;
}

void set_callback(ABM *abm, int (*callback)(double*, double[], void*),
                  double *callback_t) {
  abm->callback_t = callback_t;
  abm->callback = callback;
}

double *get_final_state(ABM *abm) {
  return abm->final_state;
}

void set_delayed_ranges(ABM *abm, int *ranges, int ranges_len) {
  free(abm->delayed_idxs);
  int idxs_len = 0;
  for (int i = 0; i < ranges_len; i += 2) {
    idxs_len += ranges[i + 1] - ranges[i];
  }
  int *idxs = (int *) malloc(idxs_len * sizeof(int));
  int k = 0;
  for (int i = 0; i < ranges_len; i += 2) {
    for (int j = ranges[i]; j < ranges[i + 1]; j++) {
      idxs[k] = j;
      k += 1;
    }
  }
  abm->delayed_idxs = idxs;
  abm->delayed_idxs_len = idxs_len;
}
