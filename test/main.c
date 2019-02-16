#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "abmd.h"
#include "plot.h"

#define DIM 4

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
  t[0] += 1 / 32.0;
//  t[0] += 1 / 16.0;
  return 1;
}

int callback_back(double *t, double *state, void *context) {
  ABMTest *abm_test = (ABMTest *) context;
  int dim = abm_test->dim;
  memcpy(&abm_test->sol_back[abm_test->i * dim], state, dim * sizeof(double));
  abm_test->i++;
  t[0] -= 1 / 32.0;
//  t[0] -= 1 / 16.0;
  return 1;
}

void calc_difference(RHSD f) {
  int order = 11;
//  double init[] = {-3844e5, 0, 0, 1023 * 3600 * 24};
  double t0 = 0;
  double t1 = 365;
  double h = 1 / 16.0;
  double delay = 0.096;
//  double delay = 0;
  int dim = DIM;
  double *init = malloc(dim * sizeof(double));
  for (int i = 0; i < dim; i += 4) {
    init[i] = -3844e5;
    init[i + 1] = 0;
    init[i + 2] = 0;
    init[i + 3] = 1023 * 3600 * 24;
  }


  int n = (int)(1 + (t1 - t0) / h);
  int sol_size = 2 * n - 1;
//  int sol_size = n;
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

  ABMD *abm = abmd_create(NULL, dim, t0, t1, h, init);
  abmd_set_f2(abm, f);
  abmd_set_delays_poly_degree(abm, 8);
  abmd_set_pointsave_poly_degree(abm, 8);
  abmd_set_delays(abm, (double[]) {delay, delay}, 2);
  abmd_set_callback(abm, callback_there, &callback_t);
  abmd_set_context(abm, &abm_test);
  abmd_set_delayed_ranges(abm, (int[]) {0, 1, 2, 3}, 4, 0);
  abmd_set_delayed_ranges(abm, (int[]) {0, 4}, 2, 1);
  abmd_set_dx_delays(abm, (int[]) {0, 1}, 2);

  abmd_run(abm);
  printf("Final: %e %e\n", abmd_get_final_state(abm)[0], abmd_get_final_state(abm)[2]);
  abmd_destroy(abm);
  printf("-----------------------------------------------------------\n");

  double *sol_reversed = malloc(sizeof(double) * sol_size * dim);

  for (int i = 0; i < sol_size; i++) {
    memcpy(&sol_reversed[i * dim], &sol[(sol_size - i - 1) * dim], sizeof(double) * dim);
    sol_reversed[i * dim + 1] *= -1;
    sol_reversed[i * dim + 3] *= -1;
  }

  abm = abmd_create(NULL, dim, t1, t0, h, &sol[(sol_size - 1) * dim]);
  abmd_set_f2(abm, f);
  abmd_set_delays_poly_degree(abm, 8);
  abmd_set_pointsave_poly_degree(abm, 8);
  abmd_set_delays(abm, (double[]) {delay, delay}, 2);
  abmd_set_callback(abm, callback_back, &callback_t);
  abmd_set_context(abm, &abm_test);
  callback_t = t1;
  abm_test.i = 0;
  abmd_set_delayed_ranges(abm, (int[]) {0, 1, 2, 3}, 4, 0);
  abmd_set_delayed_ranges(abm, (int[]) {0, 4}, 2, 1);
  abmd_set_dx_delays(abm, (int[]) {0, 1}, 2);

  abmd_run(abm);
  abmd_destroy(abm);
  double *diff = malloc(sizeof(double) * sol_size * 2);

  for (int i = 0; i < sol_size; i++) {
    double x1 = sol_reversed[i * dim];
    double x2 = sol_back[i * dim];
    double y1 = sol_reversed[i * dim + 2];
    double y2 = sol_back[i * dim + 2];
    diff[i * 2] = t0 + i * h;
    diff[i * 2 + 1] = (double) sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
  }

  plot("diff.png", diff, 2, sol_size);
  free(sol);
  free(sol_reversed);
  free(sol_back);
  free(diff);
}

void orbit(DOUBLE x[], DOUBLE xs_delayed[], DOUBLE dxs_delayed[],
           double t, DOUBLE *out, void *context) {
  int dim = 4;
  const DOUBLE G = 0.49821740236800005;
  const double m1 = 5.972e24;
  const double m2 = 7.34767309e22;
  const DOUBLE mu1 = G * m1;
  const DOUBLE mu2 = G * m2;
  const DOUBLE q = mu1 + mu2;
  DOUBLE x_delayed = xs_delayed[0];
  DOUBLE y_delayed = xs_delayed[1];

  assert(xs_delayed[0] == xs_delayed[0]);
  assert(xs_delayed[0] == xs_delayed[2]);
  assert(xs_delayed[1] == xs_delayed[4]);

  if (dxs_delayed != NULL) {
    assert(dxs_delayed[0] == dxs_delayed[0]);
    assert(dxs_delayed[0] == dxs_delayed[2]);
    assert(dxs_delayed[1] == dxs_delayed[4]);
  }

  for (int i = 0; i < DIM; i += 4) {
    DOUBLE x_ = x[i];
    DOUBLE vx = x[i + 1];
    DOUBLE y = x[i + 2];
    DOUBLE vy = x[i + 3];
    DOUBLE r = sqrt(x_ * x_ + y * y);
    out[i] = vx;
    out[i + 1] = -q * x_ / (r * r * r);
    out[i + 2] = vy;
    out[i + 3] = -q * y / (r * r * r);
    const int ae = 6371000;
    const double k2 = 0.335;
    DOUBLE c =
            -(3 * k2 * mu2 / (r * r * r)) * (1 + mu2 / mu1) * (pow(ae / r, 5));
    out[i + 1] += c * x_delayed;
    out[i + 3] += c * y_delayed;
  }
}


int main() {
  clock_t t;
  t = clock();
  calc_difference(orbit);
  t = clock() - t;
  double secs = ((double) t) / CLOCKS_PER_SEC;
//  printf("Finished in %f seconds\n", secs);
  printf("%d, %f\n", DIM, secs);
  return 0;
}
