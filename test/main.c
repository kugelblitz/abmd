#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "abm.h"
#include "plot.h"

typedef struct {
  double *callback_t;
  int i;
  int dim;
  double *sol;
  double *sol_back;
} ABMTest;

int callback_there(double *t, double *state, void *context) {
  ABMTest *abm_test = (ABMTest *) context;
  int dim = abm_test->dim;
  memcpy(&abm_test->sol[abm_test->i * dim], state, dim * sizeof(double));
//  printf("%e %e %e\n", *t, state[0], state[2]);
  abm_test->i++;
  t[0] += 1 / 32.0;
  return 1;
}

int callback_back(double *t, double *state, void *context) {
  ABMTest *abm_test = (ABMTest *) context;
  int dim = abm_test->dim;
  memcpy(&abm_test->sol_back[abm_test->i * dim], state, dim * sizeof(double));
//  printf("%e %e %e\n", *t, state[0], state[2]);
  abm_test->i++;
  t[0] -= 1 / 32.0;
  return 1;
}

void calc_difference(RHS f) {
  int order = 11;
  double init[] = {-3844e5, 0, 0, 1023 * 3600 * 24};
  double t0 = 0;
  double t1 = 5;
  double h = 1 / 16.0;
  double delay = 0;
  int dim = 4;

  int n = (int)(1 + (t1 - t0) / h);
  int sol_size = 2 * n - 1;
  double *sol = (double *) malloc(sizeof(double) * sol_size * dim);
  double *sol_back = (double *) malloc(sizeof(double) * sol_size * dim);
  double callback_t = 0;

  ABMTest abm_test = (ABMTest){
          .callback_t=&callback_t,
          .i=0,
          .dim=dim,
          .sol=sol,
          .sol_back=sol_back
  };

  ABM *abm = create_abm(f, dim, t0, t1, h, init);
  set_delays(abm, (double[]){0, delay}, 2);
  set_callback(abm, callback_there, &callback_t);
  set_context(abm, &abm_test);

  run_abm(abm);
  printf("Final: %e %e\n", get_final_state(abm)[0], get_final_state(abm)[2]);
  destroy_abm(abm);
  printf("-----------------------------------------------------------\n");

  double *sol_reversed = malloc(sizeof(double) * sol_size * dim);

  for (int i = 0; i < sol_size; i++) {
    memcpy(&sol_reversed[i * dim], &sol[(sol_size - i - 1) * dim], sizeof(double) * dim);
    sol_reversed[i * dim + 1] *= -1;
    sol_reversed[i * dim + 3] *= -1;
  }

  abm = create_abm(f, dim, t1, t0, h, &sol[(sol_size - 1) * dim]);
  set_delays(abm, (double[]){0, delay}, 2);
  set_callback(abm, callback_back, &callback_t);
  set_context(abm, &abm_test);
  callback_t = t1;
  abm_test.i = 0;

  run_abm(abm);
  destroy_abm(abm);
  double *diff = malloc(sizeof(double) * sol_size * 2);

  for (int i = 0; i < sol_size; i++) {
    double x1 = sol_reversed[i * dim];
    double x2 = sol_back[i * dim];
    double y1 = sol_reversed[i * dim + 2];
    double y2 = sol_back[i * dim + 2];
//    printf("(%.16f, %.16f) (%.16f, %.16f)\n", x1, y1, x2, y2);
    diff[i * 2] = t0 + i * h;
    diff[i * 2 + 1] = (double) sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
//    printf("%e\n", diff[i * 2 + 1]);
  }

  plot("diff.png", diff, 2, sol_size);
  free(sol);
  free(sol_reversed);
  free(sol_back);
  free(diff);
}

void orbit(DOUBLE states[], DOUBLE dotstates[], double t,
           DOUBLE *out, void *context) {
  int dim = 4;
  const DOUBLE G = 0.49821740236800005;
  const double m1 = 5.972e24;
  const double m2 = 7.34767309e22;
  const DOUBLE mu1 = G * m1;
  const DOUBLE mu2 = G * m2;
  const DOUBLE q = mu1 + mu2;
  DOUBLE x = states[0];
  DOUBLE vx = states[1];
  DOUBLE y = states[2];
  DOUBLE vy = states[3];
  DOUBLE r = sqrt(x * x + y * y);
  out[0] = vx;
  out[1] = -q * x / (r * r * r);
  out[2] = vy;
  out[3] = -q * y / (r * r * r);
  DOUBLE x_delayed = states[dim];
  DOUBLE y_delayed = states[dim + 2];
  const int ae = 6371000;
  const double k2 = 0.335;
  DOUBLE c = -(3 * k2 * mu2 / (r * r * r)) * (1 + mu2 / mu1) * (pow(ae / r, 5));
  out[1] += c * x_delayed;
  out[3] += c * y_delayed;
}



int main() {
  calc_difference(orbit);
  return 0;
}
