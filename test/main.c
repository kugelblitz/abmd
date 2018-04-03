#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
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
  abm_test->i++;
  t[0] += 2000;
  return 1;
}

int callback_back(double *t, double *state, void *context) {
  ABMTest *abm_test = (ABMTest *) context;
  int dim = abm_test->dim;
  memcpy(&abm_test->sol_back[abm_test->i * dim], state, dim * sizeof(double));
  abm_test->i++;
  t[0] -= 2000;
  return 1;
}

void calc_difference(void (*f)(DOUBLE *, DOUBLE *, double, DOUBLE *, void *)) {
  int order = 11;
  double init[] = {-3844e5, 0, 0, 1023};
  double t0 = 0;
  double t1 = 3e7;
  double h = 2000;
  double delay = 3000;
  int dim = 4;

  int n = (int)(1 + (t1 - t0) / h);
  double *sol = (double *) malloc(sizeof(double) * n * dim);
  double *sol_back = (double *) malloc(sizeof(double) * n * dim);
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
  destroy_abm(abm);
  double *sol_reversed = malloc(sizeof(double) * n * dim);

  for (int i = 0; i < n; i++) {
    memcpy(&sol_reversed[i * dim], &sol[(n - i - 1) * dim], sizeof(double) * dim);
    sol_reversed[i * dim + 1] *= -1;
    sol_reversed[i * dim + 3] *= -1;
  }

  abm = create_abm(f, dim, t1, t0, h, &sol[(n - 1) * dim]);
  set_delays(abm, (double[]){0, delay}, 2);
  set_callback(abm, callback_back, &callback_t);
  set_context(abm, &abm_test);
  callback_t = t1;
  abm_test.i = 0;

  run_abm(abm);
  destroy_abm(abm);
  double *diff = malloc(sizeof(double) * 2 * n);

  for (int i = 0; i < n; i++) {
    double x1 = sol_reversed[i * dim];
    double x2 = sol_back[i * dim];
    double y1 = sol_reversed[i * dim + 2];
    double y2 = sol_back[i * dim + 2];
    diff[i * 2] = t0 + i * h;
    diff[i * 2 + 1] = (double) sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
  }

  plot("diff.png", diff, 2, n);
  free(sol);
  free(sol_reversed);
  free(sol_back);
  free(diff);
}

void orbit(DOUBLE states[], DOUBLE dotstates[], double t,
           DOUBLE *out, void *context) {
  int dim = 4;
  const DOUBLE G = 6.67408e-11;
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