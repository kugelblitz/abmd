#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "abm.h"

DOUBLE lagrange_2d(double x, double xs[], DOUBLE ys[], int n) {
  DOUBLE y = 0;
  for (int i = 0; i < n; i++) {
    DOUBLE s = 1;
    DOUBLE t = 1;
    for (int j = 0; j < n; j++) {
      if (i == j) continue;
      s = s * (x - xs[j]);
      t = t * (xs[i] - xs[j]);
    }
    y += s / t * ys[i];
  }
	return y;
}

void lagrange(double x, double *xs, DOUBLE *ys, int dim, int n, DOUBLE *out) {
  DOUBLE *ysi = malloc(sizeof(DOUBLE) * n);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < n; j++) {
      ysi[j] = ys[j * dim + i];
    }
    out[i] = lagrange_2d(x, xs, ysi, n);
  }
  free(ysi);
}

#ifdef DEBUG

void divide(DOUBLE *qs, DOUBLE xk, int n, DOUBLE *out) {
  for (int i = 0; i < n; i++) {
    out[i] = 0;
  }
  out[n - 2] = qs[n - 1];
  for (int i = 2; i < n; i++) {
    out[n - i - 1] = qs[n - i] + xk * out[n - i];
  }
}

void compute_q(DOUBLE *xs, int n, DOUBLE *out) {
  if (n == 1) {
    out[0] = -xs[0];
    out[1] = 1;
    return;
  }
  DOUBLE *qdash = malloc(sizeof(DOUBLE) * n);
  compute_q(xs, n - 1, qdash);
  out[0] = -qdash[0] * xs[n - 1];
  out[n] = qdash[n - 1];
  for (int i = 1; i < n; i++) {
    out[i] = qdash[i - 1] - qdash[i] * xs[n - 1];
  }
  free(qdash);
}

DOUBLE evaluate_2d(DOUBLE *cs, DOUBLE x, int n) {
  DOUBLE res = 0;
  for (int i = 0; i < n; i++) {
    if (isnan(cs[i])) {
      printf("C is NaN \n");
      assert(1 == 0);
    }
    if (isnan(pow(x, i))) {
      printf("pow is NaN \n");
      assert(1 == 0);
    }
    res += cs[i] * pow(x, i);
  }
  return res;
}

void lagrange_cs_2d(DOUBLE xs[], DOUBLE ys[], int n, DOUBLE *out) {
  DOUBLE *data = malloc(sizeof(DOUBLE) * 2 * (n + 1));
  DOUBLE *q = data;
  DOUBLE *qdash = &data[n + 1];

  compute_q(xs, n, q);
  for (int i = 0; i < n; i++) {
    out[i] = 0;
  }
  for (int k = 0; k < n; k++) {
    divide(q, xs[k], n + 1, qdash);
    DOUBLE zk = evaluate_2d(qdash, xs[k], n);
    for (int i = 0; i < n; i++) {
      out[i] += ys[k] * qdash[i] / zk;
    }
  }
  free(data);
}

void lagrange_cs(DOUBLE *xs, DOUBLE *ys, int dim, int n, DOUBLE *out) {
  DOUBLE *ysi = malloc(sizeof(DOUBLE) * n);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < n; j++) {
      ysi[j] = ys[j * dim + i];
    }
    lagrange_cs_2d(xs, ysi, n, &out[i * n]);
  }
  free(ysi);
}

void evaluate(DOUBLE *cs, DOUBLE x, int dim, int n, DOUBLE *out) {
  for (int i = 0; i < dim; i++) {
    DOUBLE y = evaluate_2d(&cs[i * n], x, n);
    out[i] = y;
  }
}

void test() {
  int n = 3;
  int dim = 2;
  DOUBLE xs[3] = {0, 1, 2};
  DOUBLE ys[6] = {1, 0, 2, 1, 3, 4};
  DOUBLE *cs = malloc(sizeof(DOUBLE) * n * dim);
  lagrange_cs(xs, ys, dim, n, cs);

  DOUBLE *y = malloc(sizeof(DOUBLE) * dim);
  evaluate(cs, 3, dim, n, y);
  for (int i = 0; i < dim; i++) {
    printf("%Lf, ", y[i]);
  }
  printf("\n");
  lagrange(3, xs, ys, dim, n, y);
  for (int i = 0; i < dim; i++) {
    printf("%Lf, ", y[i]);
  }
  free(y);
  free(cs);
}

#endif