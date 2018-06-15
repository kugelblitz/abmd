#include <stdlib.h>
#include <string.h>

#include "abm.h"
#include "abm_struct.h"

#define ABM_ORDER 11

ABM *create_abm(RHS f, int dim, double t0, double t1, double h, double *init) {
  ABM *abm = (ABM *) malloc(sizeof(ABM));
  double *final_state = (double *) malloc(sizeof(double) * dim);
  char *error = (char *) malloc(256 * sizeof(char));
  *abm = (ABM) {
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
          .error=error
  };
  return abm;
}

void destroy_abm(ABM *abm) {
  free(abm->delays);
  free(abm->final_state);
  for (int i = 0; i < abm->ndelays; i++) {
    free(abm->delayed_idxs[i]);
  }
  free(abm->delayed_idxs);
  free(abm->delayed_idxs_lens);
  free(abm->error);
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

    int **idxs = (int **) malloc(sizeof(int *) * ndelays);
    int *idxs_lens = (int *) malloc(sizeof(int) * ndelays);
    memset(idxs, NULL, sizeof(int *) * ndelays);
    memset(idxs_lens, abm->dim, sizeof(int) * ndelays);
    abm->delayed_idxs = idxs;
    abm->delayed_idxs_lens = idxs_lens;
  }
}

void set_delays_poly_degree(ABM *abm, int deg) {
  abm->delays_poly_degree = deg;
}

void set_pointsave_poly_degree(ABM *abm, int deg) {
  abm->pointsave_poly_degree = deg;
}

void set_f2(ABM *abm, RHSD f2) {
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

void set_delayed_ranges(ABM *abm, int *ranges, int ranges_len, int delay_idx) {
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

char *get_last_error(ABM *abm) {
  return abm->error;
}
