#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "abmd.h"
#include "abm_struct.h"
#include "rk.h"
#include "queue.h"

static DOUBLE PREDICTOR_COEFFS[19] = {
  1.0L,
  1.0 / 2.0L,
  5.0 / 12.0L,
  3.0 / 8.0L,
  251.0 / 720.0L,
  95.0 / 288.0L,
  19087.0 / 60480.0L,
  5257.0 / 17280.0L,
  1070017.0 / 3628800.0L,
  25713.0 / 89600.0L,
  26842253.0 / 95800320.0L,
  4777223.0 / 17418240.0L,
  703604254357.0 / 2615348736000.0L,
  106364763817.0 / 402361344000.0L,
  1166309819657.0 / 4483454976000.0L,
  2.5221445e7 / 9.8402304e7L,
  8.092989203533249e15 / 3.201186852864e16L,
  8.5455477715379e13 / 3.4237292544e14L,
  1.2600467236042756559e19 / 5.109094217170944e19L
};


typedef struct {
  ABMD input;
  double rk4_h;
  DOUBLE *temp;
  Queue *queue;
  DOUBLE *xs_delayed;
  DOUBLE *xs_delayed_inner;
  DOUBLE *xs_delayed_tmp;
  DOUBLE *dxs_delayed;
  DOUBLE *rk_memory;
  DOUBLE *inner_rk_memory;
} ABMData;

void destroy_abm_data(ABMData abm_data) {
  free(abm_data.temp);
  free(abm_data.xs_delayed);
  free(abm_data.xs_delayed_inner);
  free(abm_data.xs_delayed_tmp);
  free(abm_data.dxs_delayed);
  free(abm_data.rk_memory);
  free(abm_data.inner_rk_memory);
  destroy_queue(abm_data.queue);
}

void predict(ABMData *abm_data) {

  int dim = abm_data->input.dim;
  int abm_order = abm_data->input.abm_order;
  double h = abm_data->input.h;
  Queue *queue = abm_data->queue;

  pop(queue);
  DOUBLE *prev = peek_right_x(queue);
  DOUBLE *out = push(queue);

  memset(out, 0, sizeof(DOUBLE) * dim);
  SETENV;

  DOUBLE *diffs = get_diffs_r(queue);
  for (int i = 0; i < dim; i++) {
    DOUBLE *c = PREDICTOR_COEFFS;
    for (int j = 0; j < abm_order; j++) {
      DOUBLE ch = *c++ * h;
      out[i] += *diffs++ * ch;
    }
    out[i] += prev[i];
  }
}

void correct(ABMData *abm_data, DOUBLE *x_predicted) {

  int dim = abm_data->input.dim;
  int abm_order = abm_data->input.abm_order;
  double h = abm_data->input.h;
  Queue *queue = abm_data->queue;
  DOUBLE *out = peek_right_x(queue);

  if (x_predicted == NULL) x_predicted = out;

  DOUBLE *diff = get_last_diff(queue);
  DOUBLE ch = h * PREDICTOR_COEFFS[abm_order];
  for (int k = 0; k < dim; k++) {
    out[k] = x_predicted[k] + ch * diff[k];
  }
}

void copy_delayed_states(DOUBLE x[], DOUBLE xs_delayed[], int ndelays,
                         int **idxs, int *idxs_lens) {
  for (int i = 0; i < ndelays; i++) {
    for (int j = 0; j < idxs_lens[i]; j++) {
      int idx = (idxs[i] == NULL) ? j : idxs[i][j];
      *xs_delayed++ = x[idx];
    }
  }
}

void rhs(DOUBLE x[], DOUBLE xs_delayed[], DOUBLE dxs_delayed[],
         double t, DOUBLE *out, void *abm_data) {
  ABMData *data = (ABMData *)abm_data;
  int dim = data->input.dim;
  DOUBLE *temp = data->temp;

  memset(out, 0, sizeof(DOUBLE) * dim);

  if (data->input.f1 == NULL) {
    data->input.f2(x, xs_delayed, dxs_delayed, t, out, data->input.context);
    SETENV;
    return;
  }
  
  data->input.f1(x, t, out, data->input.context);

  if (data->input.f2 == NULL) return;

  DOUBLE *out2 = &temp[0];
  memset(out2, 0, sizeof(DOUBLE) * dim);
  data->input.f2(x, xs_delayed, dxs_delayed, t, out2, data->input.context);

  SETENV;

  for (int i = 0; i < dim; i++) {
    out[i] += out2[i];
  }
}

void rhs_rk_inner(DOUBLE x[], double t, DOUBLE *out, void *abm_data) {
  ABMData *data = (ABMData *) abm_data;
  int ndelays = data->input.ndelays;
  int **idxs = data->input.delayed_idxs;
  int *idxs_lens = data->input.delayed_idxs_lens;
  DOUBLE *xs_delayed = data->xs_delayed_inner;
  copy_delayed_states(x, xs_delayed, ndelays, idxs, idxs_lens);
  rhs(x, xs_delayed, NULL, t, out, data);
}

void rhs_rk4(DOUBLE x[], double t, DOUBLE *out, void *abm_data) {

  ABMData *data = (ABMData *) abm_data;
  int dim = data->input.dim;
  int ndelays = data->input.ndelays;
  int **idxs = data->input.delayed_idxs;
  int *idxs_lens = data->input.delayed_idxs_lens;

  DOUBLE *xs_delayed = data->xs_delayed;
  DOUBLE *xs_delayed_tmp = data->xs_delayed_tmp;


  for (int i = 0; i < ndelays; i++) {
    double delay = data->input.delays[i];
    rk_step(rhs_rk_inner, -delay, t, x, dim, data, &xs_delayed_tmp[i * dim],
            NULL, &data->inner_rk_memory, METHOD_RK4);
  }

  int k = 0;
  for (int i = 0; i < ndelays; i++) {
    for (int j = 0; j < idxs_lens[i]; j++) {
      int idx = (idxs[i] == NULL) ? j : idxs[i][j];
      xs_delayed[k] = xs_delayed_tmp[i * dim + idx];
      k++;
    }
  }

  rhs(x, xs_delayed, NULL, t, out, data);
}

void get_delayed_states(ABMData *abm_data, double ti, int last_dx_known,
                        DOUBLE *x_out, DOUBLE *dx_out) {

  Queue *q = abm_data->queue;
  int ndelays = abm_data->input.ndelays;
  double *delays = abm_data->input.delays;
  int **idxs = abm_data->input.delayed_idxs;
  int *idxs_lens = abm_data->input.delayed_idxs_lens;
  int x_start_idx = 0;

  for (int i = 0; i < ndelays; i++) {
    double delay = delays[i];
    double t = ti - delay;

    evaluate_x_idxs(q, t, idxs[i], idxs_lens[i], &x_out[x_start_idx]);
    x_start_idx += idxs_lens[i];
  }

  if (dx_out == NULL) {
    return;
  }

  int dx_start_idx = 0;
  int *dx_delay_idxs = abm_data->input.dx_delays_idxs;

  for (int i = 0; i < abm_data->input.dx_delays_len; i++) {
    int delay_idx = dx_delay_idxs == NULL ? i : dx_delay_idxs[i];
    double t = ti - delays[delay_idx];
    evaluate_dx(q, t, idxs[i], idxs_lens[i], last_dx_known,
                &dx_out[dx_start_idx]);
    dx_start_idx += idxs_lens[i];
  }
}


int abmd_run(ABMD *abm) {

  int abm_order = abm->abm_order;

  if (!(1 <= abm_order && abm_order <= 19)) {
    abm->error = "ABM order must be not less than 1 and not greater than 19";
    return 1;
  }

  int delays_degree = abm->delays_poly_degree;
  int pointsave_degree = abm->pointsave_poly_degree;
  if (!(1 <= delays_degree && delays_degree <= abm_order) ||
      !(1 <= pointsave_degree && pointsave_degree <= abm_order)) {
    abm->error = "Interpolation degrees must be not less than 1 and"
                 "not greater than ABM order";
    return 1;
  }

  if (abm->f1 == NULL && abm->f2 == NULL) {
    abm->error = "Both RHSs are NULL";
    return 1;
  }

  int dim = abm->dim;
  double t0 = abm->t0;
  double t1 = abm->t1;
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
  DOUBLE *xs_delayed_tmp = (DOUBLE *) malloc(sizeof(DOUBLE) * ndelays * dim);

  int x_total_delays_len = 0;
  int dx_total_delays_len = 0;
  for (int i = 0; i < ndelays; i++) {
    x_total_delays_len += abm->delayed_idxs_lens[i];
  }
  for (int i = 0; i < abm->dx_delays_len; i++) {
    int *idxs = abm->dx_delays_idxs;
    int delay_idx = idxs == NULL ? i : idxs[i];
    dx_total_delays_len += abm->delayed_idxs_lens[delay_idx];
  }
  size_t x_delayed_size = x_total_delays_len * sizeof(DOUBLE);

  DOUBLE *xs_delayed = (DOUBLE *) malloc(x_delayed_size);
  DOUBLE *xs_delayed_inner = (DOUBLE *) malloc(x_delayed_size);
  DOUBLE *dxs_delayed = (DOUBLE *) malloc(dx_total_delays_len * sizeof(DOUBLE));

  int queue_size = abm_order + 1;
  Queue *queue = create_queue(queue_size, dim);
  qset_delays_poly_degree(queue, delays_degree);
  qset_pointsave_poly_degree(queue, pointsave_degree);
  set_t0(queue, t0);
  set_step(queue, h);

  ABMData abm_data = (ABMData){
          .input=*abm,
          .rk4_h=rk4_h,
          .temp=rhs_temp,
          .queue=queue,
          .xs_delayed=xs_delayed,
          .xs_delayed_inner=xs_delayed_inner,
          .xs_delayed_tmp=xs_delayed_tmp,
          .dxs_delayed=dxs_delayed,
          .rk_memory=NULL,
          .inner_rk_memory=NULL
  };

#ifdef __MINGW32__
  fenv_t fenv;
  fegetenv(&fenv);
#endif
  
  SETENV;

  // Setting initial conditions for RK4 solution
  for (int i = 0; i < dim; i++) {
    rk4_sol[i] = init[i];
  }

  // Doing rk4_n RK4 steps
  for (int i = 1; i < rk4_n + 1; i++) {
    double t = t0 + rk4_h * (i - 1);
    rk_step(rhs_rk4, rk4_h, t, &rk4_sol[(i - 1) * dim], dim, &abm_data,
            &rk4_sol[i * dim], &rk4_rhss[(i - 1) * dim], &abm_data.rk_memory,
            METHOD_DOPRI8);
  }

  // Computing last RHS
  int last = rk4_n * dim;
  rhs_rk4(&rk4_sol[last], t0 + rk4_h * rk4_n, &rk4_rhss[last], &abm_data);

  // Writing data from RK4 to the queue
  int k = 0;
  for (int i = 0; i < rk4_n + 1; i += RK_STEPS_IN_ABM) {
    DOUBLE *sol_address = push(queue);
    memcpy(sol_address, &rk4_sol[i * dim], dim * sizeof(DOUBLE));
    memcpy(&sol_address[dim], &rk4_rhss[i * dim], dim * sizeof(DOUBLE));
    SETENV;
    update_diffs(queue);
    swap_diffs(queue);
    k++;
  }
  free(rk4_rhss);

  int run_callback = abm->callback && abm->callback_t;
  double *callback_state = (double *) malloc(sizeof(double) * dim);
  DOUBLE *callback_state_l = (DOUBLE *) malloc(sizeof(DOUBLE) * dim);
  double *callback_t = abm->callback_t;

  SETENV;

  while (run_callback && *callback_t * hsgn <= rk4_t1 * hsgn) {
    evaluate_x_all(queue, *callback_t, callback_state_l);
    for (int i = 0; i < dim; i++) {
      callback_state[i] = (double) callback_state_l[i];
    }
    run_callback = abm->callback(callback_t, callback_state, abm->context);
    SETENV;
  }
  
  free(rk4_sol);
  free(init);

  SETENV;

  // Main ABM loop
  int start_index = rk4_i1 + 1;
  for (int i = start_index; i < n; i++) {
    double t = t0 + i * h;

    predict(&abm_data);
    DOUBLE *rhs_out = peek_right_dx(queue);
    get_delayed_states(&abm_data, t, 0, xs_delayed, dxs_delayed);
    rhs(peek_right_x(queue), xs_delayed, dxs_delayed, t, rhs_out, &abm_data);

    update_diffs(queue);
    DOUBLE *predicted = peek_right_x(queue);
    backup_last_x(queue);
    correct(&abm_data, predicted);

    get_delayed_states(&abm_data, t, 1, xs_delayed, dxs_delayed);
    rhs(peek_right_x(queue), xs_delayed, dxs_delayed, t, rhs_out, &abm_data);

    update_diffs(queue);
    restore_last_x(queue);
    correct(&abm_data, NULL);

    swap_diffs(queue);

    while (run_callback && (t - h) * hsgn < *callback_t * hsgn
                        && *callback_t * hsgn <= t * hsgn) {
      evaluate_x_all(queue, *callback_t, callback_state_l);
      for (int j = 0; j < dim; j++) {
        callback_state[j] = (double) callback_state_l[j];
      }
      run_callback = abm->callback(callback_t, callback_state, abm->context);
      SETENV;
    }
  }
  
  for (int i = 0; i < dim; i++) {
    abm->final_state[i] = (double) peek_right_x(queue)[i];
  }

  destroy_abm_data(abm_data);
  free(callback_state);
  free(callback_state_l);

  #ifdef RESTORE_FENV
  fesetenv(&fenv);
  #endif

  return 0;
}
