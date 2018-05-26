#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "abm.h"
#include "abm_struct.h"
#include "rk.h"
#include "coeffs.h"
#include "queue.h"


typedef struct {
  ABM *input;
  double rk4_h;
  DOUBLE *temp;
  Queue *queue;
  DOUBLE *states;
  DOUBLE *states_tmp;
  DOUBLE *rk_memory;
  DOUBLE *inner_rk_memory;
} ABMData;

void destroy_abm_data(ABMData abm_data) {
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
  Queue *queue = abm_data->queue;

  pop(queue);
  DOUBLE *prev = peek_right(queue);
  DOUBLE *out = push(queue);

  memset(out, 0, sizeof(DOUBLE) * dim);

  for (int j = 0; j < abm_order; j++) {
    DOUBLE *diff = get_diff(queue, j);
    DOUBLE ch = PREDICTOR_COEFFS[j] * h;
    for (int k = 0; k < dim; k++) {
      out[k] += ch * diff[k];
    }
  }

  for (int i = 0; i < dim; i++) {
    out[i] += prev[i];
  }
}

void correct(ABMData *abm_data, DOUBLE *diffs) {

  int dim = abm_data->input->dim;
  int abm_order = abm_data->input->abm_order;
  double h = abm_data->input->h;
  Queue *queue = abm_data->queue;
  DOUBLE *out = peek_right(queue);

  DOUBLE *last_diff = &diffs[abm_order * dim];
  DOUBLE ch = h * PREDICTOR_COEFFS[abm_order];
  for (int k = 0; k < dim; k++) {
    out[k] += ch * last_diff[k];
  }
}


void rhs(DOUBLE states[], DOUBLE dotstates[], double t,
         DOUBLE *out, void *abm_data) {
  ABMData *data = (ABMData *)abm_data;
  int dim = data->input->dim;
  DOUBLE *temp = data->temp;

  memset(out, 0, sizeof(DOUBLE) * dim);
  
  if (data->input->f1 != NULL) {
    data->input->f1(states, t, out, data->input->context);
    
    if (data->input->f2 != NULL) {
      DOUBLE *out2 = &temp[0];
      
      memset(out2, 0, sizeof(DOUBLE) * dim);
      data->input->f2(states, dotstates, t, out2, data->input->context);
      for (int i = 0; i < dim; i++) {
        out[i] += out2[i];
      }
    }
  }
  else if (data->input->f2 != NULL) {
    data->input->f2(states, dotstates, t, out, data->input->context);
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

void get_delayed_states(ABMData *abm_data, double ti, DOUBLE *out) {
  int dim = abm_data->input->dim;
  int ndelays = abm_data->input->ndelays;
  double *delays = abm_data->input->delays;

  for (int j = 0; j < ndelays; j++) {
    double delay = delays[j];
    double t_delayed = ti - delay;
    int all_idxs = (j == 0);

    int *idxs = all_idxs ? NULL : abm_data->input->delayed_idxs;
    int idxs_len = all_idxs ? dim : abm_data->input->delayed_idxs_len;

    evaluate_x_idxs(abm_data->queue, t_delayed, idxs, idxs_len, &out[j * dim]);
  }
}

void get_delayed_dotstates(ABMData *abm_data, double ti,
                           int last_known, DOUBLE *out) {

  if (out == NULL) return;

  int dim = abm_data->input->dim;
  int ndelays = abm_data->input->ndelays;
  double *delays = abm_data->input->delays;

  int i = 0;
  for (int j = 0; j < ndelays; j++) {
    double delay = delays[j];
    if (delay == 0) continue;
    double t_delayed = ti - delay;
    evaluate_xdot(abm_data->queue, t_delayed, abm_data->input->delayed_idxs,
                  abm_data->input->delayed_idxs_len, last_known, &out[i * dim]);
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

  int RK_STEPS_IN_ABM = 8;

  int non_zero_delays = 0;
  for (int i = 0; i < ndelays; i++) {
    if (delays[i] != 0) non_zero_delays += 1;
  }

  int n = (int)(1 + (t1 - t0) / h);
  int rk4_i1 = abm_order - 1;
  rk4_i1 += 1;  // One more step to fill the whole queue with RK data
  double rk4_h = h / (double) RK_STEPS_IN_ABM;
  int rk4_n = rk4_i1 * RK_STEPS_IN_ABM;
  double rk4_t1 = t0 + rk4_i1 * h;

  int rk_size = rk4_n + 1;  // One more block to store the initial condition
  DOUBLE *rk4_sol = (DOUBLE *) malloc(sizeof(DOUBLE) * rk_size * dim);
  DOUBLE *rk4_rhss = (DOUBLE *) malloc(sizeof(DOUBLE) * rk_size * dim);
  DOUBLE *rhs_temp = (DOUBLE *) malloc(sizeof(DOUBLE) * 2 * dim);
  DOUBLE *states_tmp = (DOUBLE *) malloc(sizeof(DOUBLE) * ndelays * dim);
  DOUBLE *states = (DOUBLE *) malloc(sizeof(DOUBLE) * (dim +
                                     abm->delayed_idxs_len * (ndelays - 1)));
  DOUBLE *dotstates = NULL;
  if (non_zero_delays > 0) {
    dotstates = (DOUBLE *) malloc(sizeof(DOUBLE) * dim *
                                  abm->delayed_idxs_len * non_zero_delays);
  }

  int queue_size = abm_order + 1;
  Queue *queue = create_queue(queue_size, 2 * dim);
  set_t0(queue, t0);
  set_step(queue, h);

  ABMData abm_data = (ABMData){
          .input=abm,
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
  for (int i = 1; i < rk4_n + 1; i++) {
    double t = t0 + rk4_h * (i - 1);
    rk_step(rhs_rk4, rk4_h, t, &rk4_sol[(i - 1) * dim], dim, 1,
            abm->delayed_idxs, abm->delayed_idxs_len, &abm_data,
            &rk4_sol[i * dim], &rk4_rhss[(i - 1) * dim], &abm_data.rk_memory,
            METHOD_DOPRI8);
  }

  // Computing last RHS
  int last = rk4_n * dim;
  rhs_rk4(&rk4_sol[last], NULL, t0 + rk4_h * (rk4_n - 1),
          &rk4_rhss[last], &abm_data);

  // Writing data from RK4 to the queue
  int k = 0;
  for (int i = 0; i < rk4_n + 1; i += RK_STEPS_IN_ABM) {
    DOUBLE *sol_address = push(queue);
    memcpy(sol_address, &rk4_sol[i * dim], dim * sizeof(DOUBLE));
    memcpy(&sol_address[dim], &rk4_rhss[i * dim], dim * sizeof(DOUBLE));
    update_diffs(queue, 0);
    k++;
  }
  free(rk4_rhss);

  int run_callback = abm->callback && abm->callback_t;
  double *callback_state = (double *) malloc(sizeof(double) * dim);
  DOUBLE *callback_state_l = (DOUBLE *) malloc(sizeof(DOUBLE) * dim);
  double *callback_t = abm->callback_t;

  while (run_callback && *callback_t * hsgn <= rk4_t1 * hsgn) {
    evaluate_x_all(queue, *callback_t, callback_state_l);
    for (int i = 0; i < dim; i++) {
      callback_state[i] = (double) callback_state_l[i];
    }
    run_callback = abm->callback(callback_t, callback_state, abm->context);
  }
  
  free(rk4_sol);
  free(init);

  DOUBLE *backup = (DOUBLE *) malloc(dim * sizeof(DOUBLE));

  // Main ABM loop
  int start_index = rk4_i1 + 1;
  for (int i = start_index; i < n; i++) {
    double t = t0 + i * h;

    predict(&abm_data);
    DOUBLE *rhs_out = &peek_right(queue)[dim];
    get_delayed_states(&abm_data, t, states);
    get_delayed_dotstates(&abm_data, t, 0, dotstates);
    rhs(states, dotstates, t, rhs_out, &abm_data);

    DOUBLE *diffs = update_diffs(queue, 1);
    memcpy(backup, peek_right(queue), dim * sizeof(DOUBLE));
    correct(&abm_data, diffs);

    get_delayed_states(&abm_data, t, states);
    get_delayed_dotstates(&abm_data, t, 1, dotstates);
    rhs(states, dotstates, t, rhs_out, &abm_data);

    diffs = update_diffs(queue, 0);
    memcpy(peek_right(queue), backup, dim * sizeof(DOUBLE));
    correct(&abm_data, diffs);

    while (run_callback && (t - h) * hsgn < *callback_t * hsgn
                        && *callback_t * hsgn <= t * hsgn) {
      evaluate_x_all(queue, *callback_t, callback_state_l);
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
  free(backup);
  free(callback_state);
  free(callback_state_l);
}
