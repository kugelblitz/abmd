#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "abm.h"
#include "abm_struct.h"
#include "rk.h"
#include "coeffs.h"
#include "poly.h"
#include "queue.h"


typedef struct {
  ABM *input;
  DOUBLE *predictor_coeffs;
  DOUBLE *corrector_coeffs;
  double *interp_xs;
  DOUBLE *interp_ys;
  double *extrap_xs;
  DOUBLE *extrap_ys;
  DOUBLE *lagrange_data;
  double rk4_h;
  DOUBLE *temp;
  Queue *queue;
  DOUBLE *states;
  DOUBLE *states_tmp;
  DOUBLE *rk_memory;
  DOUBLE *inner_rk_memory;
} ABMData;

void destroy_abm_data(ABMData abm_data) {
  free(abm_data.predictor_coeffs);
  free(abm_data.corrector_coeffs);
  free(abm_data.interp_xs);
  free(abm_data.interp_ys);
  free(abm_data.extrap_xs);
  free(abm_data.extrap_ys);
  free(abm_data.lagrange_data);
  free(abm_data.temp);
  free(abm_data.states);
  free(abm_data.states_tmp);
  free(abm_data.rk_memory);
  free(abm_data.inner_rk_memory);
  destroy_queue(abm_data.queue);
}

void predict(ABMData *abm_data) {

  int dim = abm_data->input->dim;
  int abm_order = abm_data->input->abm_order;
  double h = abm_data->input->h;
  DOUBLE *coeffs = abm_data->predictor_coeffs;
  Queue *queue = abm_data->queue;
  int q_size = get_capacity(queue);
  DOUBLE *temp = abm_data->temp;

  pop(queue);
  DOUBLE *prev = peek_right(queue);
  DOUBLE *out = push(queue);
  memcpy(out, prev, sizeof(DOUBLE) * dim);

  memset(temp, 0, sizeof(DOUBLE) * dim);

  for (int j = 0; j < abm_order; j++) {
    DOUBLE *rhs = &get(queue, q_size - j - 2)[dim];
    DOUBLE c = coeffs[j];
    for (int k = 0; k < dim; k++) {
      temp[k] += c * h * rhs[k];
    }
  }

  for (int i = 0; i < dim; i++) {
    out[i] += temp[i];
  }
}

void correct(ABMData *abm_data) {

  int dim = abm_data->input->dim;
  int abm_order = abm_data->input->abm_order;
  double h = abm_data->input->h;
  DOUBLE *coeffs = abm_data->corrector_coeffs;
  Queue *queue = abm_data->queue;
  int q_size = get_capacity(queue);
  DOUBLE *prev = get(queue, q_size - 2);
  DOUBLE *out = peek_right(queue);
  DOUBLE *temp = abm_data->temp;

  memcpy(out, prev, sizeof(DOUBLE) * dim);

  memset(temp, 0, sizeof(DOUBLE) * dim);

  for (int j = 0; j < abm_order; j++) {
    DOUBLE *rhs = &get(queue, q_size - j - 1)[dim];
    DOUBLE c = coeffs[j];
    for (int k = 0; k < dim; k++) {
      temp[k] += c * h * rhs[k];
    }
  }

  for (int i = 0; i < dim; i++) {
    out[i] += temp[i];
  }
}


void rhs(DOUBLE states[], DOUBLE dotstates[], double t,
         DOUBLE *out, void *abm_data) {
  ABMData *data = (ABMData *)abm_data;
  int dim = data->input->dim;
  DOUBLE *temp = data->temp;

  memset(temp, 0, sizeof(DOUBLE) * 2 * dim);

  DOUBLE *out1 = &temp[0];
  DOUBLE *out2 = &temp[dim];

  if (data->input->f1 != NULL) {
    data->input->f1(states, t, out1, data->input->context);
  }
  if (data->input->f2 != NULL) {
    data->input->f2(states, dotstates, t, out2, data->input->context);
  }
  for (int i = 0; i < dim; i++) {
    out[i] = out1[i] + out2[i];
  }
}

void rhs_rk4(DOUBLE state[], DOUBLE dotstates[], double t,
             DOUBLE *out, void *abm_data) {

  ABMData *data = (ABMData *)abm_data;
  int dim = data->input->dim;
  int ndelays = data->input->ndelays;
  int *delayed_idxs = data->input->delayed_idxs;
  int delayed_idxs_len = data->input->delayed_idxs_len;

  DOUBLE *states_tmp = data->states_tmp;
  DOUBLE *states = data->states;

  for (int i = 0; i < ndelays; i++) {
    double delay = data->input->delays[i];
    if (delay == 0) {
      memcpy(&states_tmp[i * dim], state, dim * sizeof(DOUBLE));
      continue;
    }
    rk_step(rhs, -delay, t, state, dim, ndelays, delayed_idxs, delayed_idxs_len,
            data, &states_tmp[i * dim], NULL, &data->inner_rk_memory, METHOD_RK4);
  }

  memcpy(states, states_tmp, dim * sizeof(DOUBLE));

  for (int i = 1; i < ndelays; i++) {
    for (int j = 0; j < delayed_idxs_len; j++) {
      int idx = delayed_idxs == NULL ? j : delayed_idxs[j];
      states[dim + (i - 1) * delayed_idxs_len + j] = states_tmp[i * dim + idx];
    }
  }

  rhs(states, NULL, t, out, data);
}

void get_state_at_time(ABMData *abm_data, double t, double t_last,
                       int derivative, int is_last_defined, int idxs_only,
                       DOUBLE *out) {
  Queue *queue = abm_data->queue;
  int q_last_idx = get_capacity(queue) - 1;
  int dim = abm_data->input->dim;
  double h = abm_data->input->h;

  int state_idx = 0;
  if (derivative == 1) {
    state_idx = dim;
    if (is_last_defined == 0) {
      t_last -= h;
      q_last_idx -= 1;
    }
  }

  double steps_diff = (t - t_last) / h;

  int *idxs = abm_data->input->delayed_idxs;
  int idxs_len = abm_data->input->delayed_idxs_len;
  int idim = (idxs_only && idxs != NULL) ? idxs_len : dim;

  if (steps_diff <= 0) {
    if (!fmod(steps_diff, 1)) {
      DOUBLE *state = get(queue, q_last_idx + (int) steps_diff);
      for (int i = 0; i < idim; i++) {
        int idx = (idxs_only && idxs != NULL) ? idxs[i] : i;
        out[i] = state[idx];
      }
      return;
    }
    // Interpolation
    int order = abm_data->input->interpolation_order;
    int points_number = order + 1;
    int right_i = q_last_idx + (int) steps_diff;
    int left_i = right_i - points_number + 1;
    double right_t = t_last + h * (int) steps_diff;
    double left_t = right_t - (points_number - 1) * h;
    if (left_i < 0) {
      left_t -= (left_i * h);
      left_i = 0;
    }
//    This was supposed to improve interpolation precision, but only seems
//        to make it worse. More over, it cannot easily be used in current
//        function, so I'll just leave it like that for now.
//    if (delay > h) {
//        left_i += 1;
//        left_t += h;
//    }
    if (abm_data->interp_xs == NULL) {
      abm_data->interp_xs = (double *) malloc(sizeof(double) * points_number);
      abm_data->interp_ys = (DOUBLE *) malloc(sizeof(DOUBLE) *
                                              points_number * dim);
    }
    double *xs = abm_data->interp_xs;
    DOUBLE *ys = abm_data->interp_ys;

    for (int ii = 0; ii < points_number; ii++) {
      xs[ii] = left_t + ii * h;
      DOUBLE *sol = &get(queue, left_i + ii)[state_idx];
      for (int jj = 0; jj < idim; jj++) {
        int idx = (idxs_only && idxs != NULL) ? idxs[jj] : jj;
        ys[ii * idim + jj] = sol[idx];
      }
    }
    lagrange(t, xs, ys, idim, points_number, &abm_data->lagrange_data, out);
    return;
  }
  // Extrapolation
  int points_number = abm_data->input->extrapolation_order + 1;
  int left_i = q_last_idx - points_number + 1;
  double left_t = t_last - (points_number - 1) * h;
  if (abm_data->extrap_xs == NULL) {
    abm_data->extrap_xs = (double *) malloc(sizeof(double) * points_number);
    abm_data->extrap_ys = (DOUBLE *) malloc(sizeof(DOUBLE) *
                                            points_number * dim);
  }
  double *xs = abm_data->extrap_xs;
  DOUBLE *ys = abm_data->extrap_ys;
  for (int ii = 0; ii < points_number; ii++) {
    xs[ii] = left_t + ii * h;
    DOUBLE *sol = &get(queue, left_i + ii)[state_idx];
    for (int jj = 0; jj < idim; jj++) {
      int idx = (idxs_only && idxs != NULL) ? idxs[jj] : jj;
      ys[ii * idim + jj] = sol[idx];
    }
  }
  lagrange(t, xs, ys, idim, points_number, &abm_data->lagrange_data, out);
}

void get_delayed_states(ABMData *abm_data, double ti, double t_last,
                        DOUBLE *out) {
  int dim = abm_data->input->dim;
  int ndelays = abm_data->input->ndelays;
  double *delays = abm_data->input->delays;

  for (int j = 0; j < ndelays; j++) {
    double delay = delays[j];
    double t_delayed = ti - delay;
    int idxs_only = j != 0;
    get_state_at_time(abm_data, t_delayed, t_last, 0, -1,
                      idxs_only, &out[j * dim]);
  }
}

void get_delayed_dotstates(ABMData *abm_data, double ti, double t_last,
                           int is_last_rhs_defined, DOUBLE *out) {

  if (out == NULL) return;

  int dim = abm_data->input->dim;
  int ndelays = abm_data->input->ndelays;
  double *delays = abm_data->input->delays;

  int i = 0;
  for (int j = 0; j < ndelays; j++) {
    double delay = delays[j];
    if (delay == 0) continue;
    double t_delayed = ti - delay;
    get_state_at_time(abm_data, t_delayed, t_last, 1, is_last_rhs_defined, 1,
                      &out[i * dim]);
    i += 1;
  }
}

void run_abm(ABM *abm) {

  int dim = abm->dim;
  int abm_order = abm->abm_order;
  double t0 = abm->t0;
  double t1 = abm->t1;
  double *delays = abm->delays;
  int ndelays = abm->ndelays;

  int hsgn = (t1 > t0) - (t1 < t0);
  abm->h *= hsgn;
  double h = abm->h;

  DOUBLE *init = (DOUBLE *) malloc(sizeof(DOUBLE) * dim);
  for (int i = 0; i < dim; i++) {
    init[i] = abm->init[i];
  }

  if (abm->init_call != NULL) {
    abm->init_call(init, abm->context);
  }

  DOUBLE *predictor_coeffs = (DOUBLE *) malloc(sizeof(DOUBLE) * abm_order);
  DOUBLE *corrector_coeffs = (DOUBLE *) malloc(sizeof(DOUBLE) * abm_order);
  get_predictor_coeffs(abm_order, predictor_coeffs);
  get_corrector_coeffs(abm_order, corrector_coeffs);

  int RK_STEPS_IN_ABM = 8;

  int interpolation_order = abm->interpolation_order;

  int non_zero_delays = 0;
  double max_positive_delay = 0;
  double max_negative_delay = 0;
  for (int i = 0; i < ndelays; i++) {
    if (delays[i] != 0) non_zero_delays += 1;
    if (delays[i] > max_positive_delay) {
      max_positive_delay = delays[i];
    }
    if (delays[i] < max_negative_delay) {
      max_negative_delay = delays[i];
    }
  }

  int n = (int)(1 + (t1 - t0) / h);
  int extra_steps = fmod(max_positive_delay, h) ?
                    (int)(max_positive_delay / h) + 1 :
                    (int)(max_positive_delay / h);
  if (max_positive_delay != 0) {
    extra_steps += interpolation_order - 1;
  }
  int rk4_i1 = abm_order - 1 + extra_steps;
  double rk4_h = h / (double) RK_STEPS_IN_ABM;
  int rk4_n = 1 + rk4_i1 * RK_STEPS_IN_ABM;
  double rk4_t1 = t0 + rk4_i1 * h;


  DOUBLE *rk4_sol = (DOUBLE *) malloc(sizeof(DOUBLE) * rk4_n * dim);
  DOUBLE *rk4_rhss = (DOUBLE *) malloc(sizeof(DOUBLE) * rk4_n * dim);
  DOUBLE *rhs_temp = (DOUBLE *) malloc(sizeof(DOUBLE) * 2 * dim);
  DOUBLE *states_tmp = (DOUBLE *) malloc(sizeof(DOUBLE) * ndelays * dim);
  DOUBLE *states = (DOUBLE *) malloc(sizeof(DOUBLE) * (dim +
                                     abm->delayed_idxs_len * (ndelays - 1)));
  DOUBLE *dotstates = NULL;
  if (non_zero_delays > 0) {
    dotstates = (DOUBLE *) malloc(sizeof(DOUBLE) * dim *
                                  abm->delayed_idxs_len * non_zero_delays);
  }

  int queue_size = (int) ceil(rk4_n / (double) RK_STEPS_IN_ABM);
  Queue *queue = create_queue(queue_size, 2 * dim);

  ABMData abm_data = (ABMData){
          .input=abm,
          .predictor_coeffs=predictor_coeffs,
          .corrector_coeffs=corrector_coeffs,
          .interp_xs=NULL,
          .interp_ys=NULL,
          .extrap_xs=NULL,
          .extrap_ys=NULL,
          .lagrange_data=NULL,
          .rk4_h=rk4_h,
          .temp=rhs_temp,
          .queue=queue,
          .states=states,
          .states_tmp=states_tmp,
          .rk_memory=NULL,
          .inner_rk_memory=NULL
  };

  // Setting initial conditions for RK4 solution
  for (int i = 0; i < dim; i++) {
    rk4_sol[i] = init[i];
  }

  // Doing rk4_n RK4 steps
  for (int i = 1; i < rk4_n; i++) {
    double t = t0 + rk4_h * (i - 1);
    rk_step(rhs_rk4, rk4_h, t, &rk4_sol[(i - 1) * dim], dim, 1,
            abm->delayed_idxs, abm->delayed_idxs_len, &abm_data,
            &rk4_sol[i * dim], &rk4_rhss[(i - 1) * dim], &abm_data.rk_memory,
            METHOD_DOPRI8);
  }

  // Computing last RHS
  int last = (rk4_n - 1) * dim;
  rhs_rk4(&rk4_sol[last], NULL, t0 + rk4_h * (rk4_n - 1),
          &rk4_rhss[last], &abm_data);

  // Writing data from RK4 to the queue
  int k = 0;
  for (int i = 0; i < rk4_n; i += RK_STEPS_IN_ABM) {
    DOUBLE *sol_address = push(queue);
    memcpy(sol_address, &rk4_sol[i * dim], dim * sizeof(DOUBLE));
    memcpy(&sol_address[dim], &rk4_rhss[i * dim], dim * sizeof(DOUBLE));
    k++;
  }
  free(rk4_rhss);

  int run_callback = abm->callback && abm->callback_t;
  double *callback_state = (double *) malloc(sizeof(double) * dim);
  DOUBLE *callback_state_l = (DOUBLE *) malloc(sizeof(DOUBLE) * dim);
  double *callback_t = abm->callback_t;

  while (run_callback && *callback_t * hsgn <= rk4_t1 * hsgn) {
    get_state_at_time(&abm_data, *callback_t, rk4_t1,
                      0, -1, 0, callback_state_l);
    for (int i = 0; i < dim; i++) {
      callback_state[i] = (double) callback_state_l[i];
    }
    run_callback = abm->callback(callback_t, callback_state, abm->context);
  }
  
  free(rk4_sol);
  free(init);


  // If ABM_ORDER = 2: rk4_t1 = 2000, t = 4000
  int start_index = rk4_i1 + 1;
  for (int i = start_index; i < n; i++) {
    double t = t0 + i * h;

    predict(&abm_data);
    DOUBLE *rhs_out = &peek_right(queue)[dim];
    get_delayed_states(&abm_data, t, t, states);
    get_delayed_dotstates(&abm_data, t, t, 0, dotstates);
    rhs(states, dotstates, t, rhs_out, &abm_data);

    correct(&abm_data);
    get_delayed_states(&abm_data, t, t, states);
    get_delayed_dotstates(&abm_data, t, t, 1, dotstates);
    rhs(states, dotstates, t, rhs_out, &abm_data);

    correct(&abm_data);

    while (run_callback && (t - h) * hsgn < *callback_t * hsgn
                        && *callback_t * hsgn <= t * hsgn) {
      get_state_at_time(&abm_data, *callback_t, t, 0, -1, 0, callback_state_l);
      for (int j = 0; j < dim; j++) {
        callback_state[j] = (double) callback_state_l[j];
      }
      run_callback = abm->callback(callback_t, callback_state, abm->context);
    }
  }
  
  for (int i = 0; i < dim; i++) {
    abm->final_state[i] = (double) peek_right(queue)[i];
  }

  destroy_abm_data(abm_data);
  free(dotstates);
  free(callback_state);
  free(callback_state_l);
}
