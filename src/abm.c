#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>

#include "abm.h"
#include "abm_struct.h"
#include "coeffs.h"
#include "poly.h"
#include "queue.h"


typedef struct {
  ABM *input;
  DOUBLE *predictor_coeffs;
  DOUBLE *corrector_coeffs;
  double rk4_h;
  DOUBLE *rhs_temp;
  Queue *queue;
} ABMData;


void predict(ABMData *abm_data) {

  int dim = abm_data->input->dim;
  int abm_order = abm_data->input->abm_order;
  double h = abm_data->input->h;
  DOUBLE *coeffs = abm_data->predictor_coeffs;
  Queue *queue = abm_data->queue;
  int q_size = get_capacity(queue);

  pop(queue);
  DOUBLE *prev = peek_right(queue);
  DOUBLE *out = push(queue);
  memcpy(out, prev, sizeof(DOUBLE) * dim);

  for (int j = 0; j < abm_order; j++) {
    DOUBLE *rhs = &get(queue, q_size - j - 2)[dim];
    DOUBLE c = coeffs[j];
    for (int k = 0; k < dim; k++) {
      out[k] += c * h * rhs[k];
    }
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

  memcpy(out, prev, sizeof(DOUBLE) * dim);


  for (int j = 0; j < abm_order; j++) {
    DOUBLE *rhs = &get(queue, q_size - j - 1)[dim];
    DOUBLE c = coeffs[j];
    for (int k = 0; k < dim; k++) {
      out[k] += c * h * rhs[k];
    }
  }
}

void multiply_vector(DOUBLE *vec, DOUBLE c, int dim) {
  for (int i = 0; i < dim; i++) {
    vec[i] = vec[i] * c;
  }
}

void rk4_step(void (*f)(double, DOUBLE[], void*, DOUBLE*), double h, double t,
              DOUBLE *state, int dim, int ndelays, void *context, DOUBLE *out) {

  DOUBLE *data = malloc(sizeof(DOUBLE) * 4 * dim * (1 + ndelays));

  DOUBLE *ks = data;
  DOUBLE *k1 = ks;
  DOUBLE *k2 = &ks[dim];
  DOUBLE *k3 = &ks[dim * 2];
  DOUBLE *k4 = &ks[dim * 3];

  DOUBLE *ins = &data[dim * 4];
  DOUBLE *in1 = ins;
  DOUBLE *in2 = &ins[dim * ndelays];
  DOUBLE *in3 = &ins[dim * ndelays * 2];
  DOUBLE *in4 = &ins[dim * ndelays * 3];

  for (int i = 0; i < ndelays; i++) {
    memcpy(&in1[i * dim], state, dim * sizeof(DOUBLE));
  }
  f(t, in1, context, k1);
  multiply_vector(k1, h, dim);

  for (int i = 0; i < ndelays; i++) {
    for (int j = 0; j < dim; j++) {
      in2[i * dim + j] = state[j] + k1[j] / 2;
    }
  }
  f(t + h / 2, in2, context, k2);
  multiply_vector(k2, h, dim);

  for (int i = 0; i < ndelays; i++) {
    for (int j = 0; j < dim; j++) {
      in3[i * dim + j] = state[j] + k2[j] / 2;
    }
  }
  f(t + h / 2, in3, context, k3);
  multiply_vector(k3, h, dim);

  for (int i = 0; i < ndelays; i++) {
    for (int j = 0; j < dim; j++) {
      in4[i * dim + j] = state[j] + k3[j];
    }
  }
  f(t + h, in4, context, k4);
  multiply_vector(k4, h, dim);

  for (int i = 0; i < dim; i++) {
    out[i] = state[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
  }
  free(data);
}

void rhs(double t, DOUBLE states[], void *abm_data, DOUBLE *out) {

  ABMData *data = (ABMData *)abm_data;
  int dim = data->input->dim;
  DOUBLE *temp = data->rhs_temp;

  memset(temp, 0, sizeof(DOUBLE) * 2 * dim);

  DOUBLE *out1 = &temp[0];
  DOUBLE *out2 = &temp[dim];

  if (data->input->f1 != NULL) {
    data->input->f1(t, states, data->input->context, out1);
  }
  if (data->input->f2 != NULL) {
    data->input->f2(t, states, data->input->context, out2);
  }
  for (int i = 0; i < dim; i++) {
    out[i] = out1[i] + out2[i];
  }
}

void rhs_rk4(double t, DOUBLE *state, void *abm_data, DOUBLE *out) {

  ABMData *data = (ABMData *)abm_data;
  int dim = data->input->dim;
  int ndelays = data->input->ndelays;
  DOUBLE *states = malloc(sizeof(DOUBLE) * ndelays * dim);
  for (int i = 0; i < ndelays; i++) {
    double delay = data->input->delays[i];
    if (delay == 0) {
      memcpy(&states[i * dim], state, dim * sizeof(DOUBLE));
      continue;
    }
    rk4_step(rhs, -delay, t, state, dim, ndelays,
         data, &states[i * dim]);
  }
  rhs(t, states, data, out);
  free(states);
}

void get_state_at_time(ABMData *abm_data, double t, double t_last, DOUBLE *out) {
  Queue *queue = abm_data->queue;
  int q_size = get_capacity(queue);
  int dim = abm_data->input->dim;
  double h = abm_data->input->h;
  double steps_diff = (t - t_last) / h;
  if (steps_diff <= 0) {
    if (!fmod(steps_diff, 1)) {
      DOUBLE *state = get(queue, q_size - 1 + (int) steps_diff);
      for (int i = 0; i < dim; i++) {
        out[i] = state[i];
      }
      return;
    }
    // Interpolation
    int order = abm_data->input->interpolation_order;
    int points_number = order + 1;
    int right_i = q_size - 1 + (int) steps_diff;
    int left_i = right_i - points_number + 1;
    double right_t = t_last + h * (int) steps_diff;
    double left_t = right_t - (points_number - 1) * h;
//        This was supposed to improve interpolation precision, but only seems
//            to make it worse. More over, it cannot easily be used in current
//            function, so I'll just leave it like that for now.
//        if (delay > h) {
//            left_i += 1;
//            left_t += h;
//        }
    double *xs = malloc(sizeof(double) * points_number);
    DOUBLE *ys = malloc(sizeof(DOUBLE) * points_number * dim);
    for (int ii = 0; ii < points_number; ii++) {
      xs[ii] = left_t + ii * h;
      DOUBLE *sol = get(queue, left_i + ii);
      for (int jj = 0; jj < dim; jj++) {
        ys[ii * dim + jj] = sol[jj];
      }
    }
    lagrange(t, xs, ys, dim, points_number, out);
    return;
  }
  // Extrapolation
  int points_number = abm_data->input->extrapolation_order + 1;
  int left_i = q_size - points_number;
  double left_t = t_last - (points_number - 1) * h;
  double *xs = malloc(sizeof(double) * points_number);
  DOUBLE *ys = malloc(sizeof(DOUBLE) * points_number * dim);
  for (int ii = 0; ii < points_number; ii++) {
    xs[ii] = left_t + ii * h;
    DOUBLE *sol = get(queue, left_i + ii);
    for (int jj = 0; jj < dim; jj++) {
      ys[ii * dim + jj] = sol[jj];
    }
  }
  lagrange(t, xs, ys, dim, points_number, out);
}

DOUBLE *get_delayed_states(double ti, double t_last, ABMData *abm_data) {

  int dim = abm_data->input->dim;
  int ndelays = abm_data->input->ndelays;
  double *delays = abm_data->input->delays;

  DOUBLE *states = malloc(sizeof(DOUBLE) * dim * ndelays);
  for (int j = 0; j < ndelays; j++) {
    double delay = delays[j];
    double t_delayed = ti - delay;
    get_state_at_time(abm_data, t_delayed, t_last, &states[j * dim]);
  }
  return states;
}

void run_abm(ABM *abm) {

  int dim = abm->dim;
  int abm_order = abm->abm_order;
  double t0 = abm->t0;
  double t1 = abm->t1;
  double h = abm->h;
  double *delays = abm->delays;
  int ndelays = abm->ndelays;

  if (t1 < t0) {
    h = -h;
    abm->h = h;
  }

  int hsgn = (h > 0) - (h < 0);

  DOUBLE *init = malloc(sizeof(DOUBLE) * dim);
  for (int i = 0; i < dim; i++) {
    init[i] = abm->init[i];
  }

  if (abm->init_call != NULL) {
    abm->init_call(init, abm->context);
  }

  DOUBLE *coeffs = (DOUBLE *) malloc(sizeof(DOUBLE) * 2 * abm_order);
  DOUBLE *predictor_coeffs = &coeffs[0];
  DOUBLE *corrector_coeffs = &coeffs[abm_order];
  get_predictor_coeffs(abm_order, predictor_coeffs);
  get_corrector_coeffs(abm_order, corrector_coeffs);

  int RK_STEPS_IN_ABM = 4;

  int interpolation_order = abm->interpolation_order;
  int extrapolation_order = abm->extrapolation_order;

  double max_positive_delay = 0;
  double max_negative_delay = 0;
  for (int i = 0; i < ndelays; i++) {
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
  DOUBLE *rhs_temp = malloc(sizeof(DOUBLE) * 2 * dim);

  int queue_size = (int) ceil(rk4_n / (double) RK_STEPS_IN_ABM);
  Queue *queue = create_queue(queue_size, 2 * dim);

  ABMData abm_data = (ABMData){
          .input=abm,
          .predictor_coeffs=predictor_coeffs,
          .corrector_coeffs=corrector_coeffs,
          .rk4_h=rk4_h,
          .rhs_temp=rhs_temp,
          .queue=queue
  };

  // Setting initial conditions for RK4 solution
  for (int i = 0; i < dim; i++) {
    rk4_sol[i] = init[i];
  }

  // Doing rk4_n RK4 steps
  for (int i = 1; i < rk4_n; i++) {
    double t = t0 + rk4_h * i;
    rk4_step(rhs_rk4, rk4_h, t, &rk4_sol[(i - 1) * dim], dim, 1,
             &abm_data, &rk4_sol[i * dim]);
  }

  // Writing data from RK4 to the queue
  int k = 0;
  for (int i = 0; i < rk4_n; i += RK_STEPS_IN_ABM) {
    DOUBLE *sol_address = push(queue);
    memcpy(sol_address, &rk4_sol[i * dim], dim * sizeof(DOUBLE));
    k++;
  }
  free(rk4_sol);

  int run_callback = abm->callback && abm->callback_t;
  double *callback_state = malloc(sizeof(double) * dim);
  DOUBLE *callback_state_l = malloc(sizeof(DOUBLE) * dim);
  double *callback_t = abm->callback_t;

  while (run_callback && *callback_t * hsgn <= rk4_t1 * hsgn) {
    get_state_at_time(&abm_data, *callback_t, rk4_t1, callback_state_l);
    for (int i = 0; i < dim; i++) {
      callback_state[i] = (double) callback_state_l[i];
    }
    run_callback = abm->callback(callback_t, callback_state, abm->context);
  }

  // Initializing right-hand sides
  for (int i = extra_steps; i < k; i++) {
    DOUBLE *states = get_delayed_states(t0 + i * h, rk4_t1, &abm_data);
    DOUBLE *address = &get(queue, i)[dim];
    rhs(t0 + i * h, states, &abm_data, address);
    free(states);
  }

  // If ABM_ORDER = 2: rk4_t1 = 2000, t = 4000
  int start_index = rk4_i1 + 1;
  for (int i = start_index; i < n; i++) {
    double t = t0 + i * h;

    predict(&abm_data);
    DOUBLE *rhs_out = &peek_right(queue)[dim];
    DOUBLE *states = get_delayed_states(t, t, &abm_data);
    rhs(t, states, &abm_data, rhs_out);

    free(states);

    correct(&abm_data);
    states = get_delayed_states(t, t, &abm_data);
    rhs(t, states, &abm_data, rhs_out);

    free(states);

    if (run_callback) {
      if ((t - h) * hsgn < *callback_t * hsgn  && *callback_t * hsgn <= t * hsgn) {
        get_state_at_time(&abm_data, *callback_t, t, callback_state_l);
        for (int j = 0; j < dim; j++) {
          callback_state[j] = (double) callback_state_l[j];
        }
        run_callback = abm->callback(callback_t, callback_state, abm->context);
      }
    }
  }

  for (int i = 0; i < dim; i++) {
    abm->final_state[i] = (double) peek_right(queue)[i];
  }

  destroy_queue(queue);
  free(coeffs);
  free(init);
  free(callback_state);
  free(callback_state_l);
}
