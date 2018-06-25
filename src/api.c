#include <stdlib.h>
#include <string.h>

#include "abmd.h"
#include "abm_struct.h"

#define ABM_ORDER 11

ABMD *abmd_create(RHS f, int dim, double t0, double t1, double h, double *init) {
  ABMD *abm = (ABMD *) malloc(sizeof(ABMD));
  double *final_state = (double *) malloc(sizeof(double) * dim);
  char *error = (char *) malloc(256 * sizeof(char));
  *abm = (ABMD) {
          .f1=f,
          .f2=NULL,
          .dim=dim,
          .t0=t0,
          .t1=t1,
          .h=h,
          .init=init,
          .delays=NULL,
          .ndelays=0,
          .abm_order=ABM_ORDER,
          .delays_poly_degree=ABM_ORDER,
          .pointsave_poly_degree=ABM_ORDER,
          .final_state=final_state,
          .context=NULL,
          .init_call=NULL,
          .callback_t=NULL,
          .callback=NULL,
          .delayed_idxs=NULL,
          .delayed_idxs_lens=NULL,
          .dx_delays_idxs=NULL,
          .dx_delays_len=0,
          .error=error
  };
  return abm;
}

void abmd_destroy(ABMD *abm) {
  free(abm->delays);
  free(abm->final_state);
  for (int i = 0; i < abm->ndelays; i++) {
    free(abm->delayed_idxs[i]);
  }
  free(abm->delayed_idxs);
  free(abm->delayed_idxs_lens);
  free(abm->dx_delays_idxs);
  free(abm->error);
  free(abm);
}

void abmd_set_order(ABMD *abm, int order) {
  abm->abm_order = order;
}

void abmd_set_delays(ABMD *abm, double *delays, int ndelays) {
  if (abm->delays != NULL) {
    free(abm->delays);
    abm->delays = NULL;
  }
  abm->ndelays = ndelays;
  abm->dx_delays_len = ndelays;
  if (ndelays > 0) {
    abm->delays = malloc(sizeof(double) * ndelays);
    memcpy(abm->delays, delays, sizeof(double) * ndelays);

    int **idxs = (int **) malloc(sizeof(int *) * ndelays);
    int *idxs_lens = (int *) malloc(sizeof(int) * ndelays);
    for (int i = 0; i < ndelays; i++) {
      idxs[i] = NULL;
      idxs_lens[i] = abm->dim;
    }
    abm->delayed_idxs = idxs;
    abm->delayed_idxs_lens = idxs_lens;
  }
}

void abmd_set_delays_poly_degree(ABMD *abm, int deg) {
  abm->delays_poly_degree = deg;
}

void abmd_set_pointsave_poly_degree(ABMD *abm, int deg) {
  abm->pointsave_poly_degree = deg;
}

void abmd_set_f2(ABMD *abm, RHSD f2) {
  abm->f2 = f2;
}

void abmd_set_context(ABMD *abm, void *context) {
  abm->context = context;
}

void abmd_set_init_call(ABMD *abm, void (*init_call)(DOUBLE[], void *)) {
  abm->init_call = init_call;
}

void abmd_set_callback(ABMD *abm, int (*callback)(double *, double[], void *),
                       double *callback_t) {
  abm->callback_t = callback_t;
  abm->callback = callback;
}

double *abmd_get_final_state(ABMD *abm) {
  return abm->final_state;
}

void abmd_set_delayed_ranges(ABMD *abm, int *ranges, int ranges_len,
                             int delay_idx) {
  free(abm->delayed_idxs[delay_idx]);
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
  abm->delayed_idxs[delay_idx] = idxs;
  abm->delayed_idxs_lens[delay_idx] = idxs_len;
}

void abmd_set_dx_delays(ABMD *abm, int *idxs, int idxs_len) {
  int *_idxs = (int *) malloc(idxs_len * sizeof(int));
  memcpy(_idxs, idxs, idxs_len * sizeof(int));
  abm->dx_delays_idxs = _idxs;
  abm->dx_delays_len = idxs_len;
}

char *abmd_get_last_error(ABMD *abm) {
  return abm->error;
}
