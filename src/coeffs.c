#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "abm.h"


DOUBLE PREDICTOR_COEFFS[19] = {
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


DOUBLE _MOULTON_COEFFS[19] = {
   1,
  -1 / 2.0L,
  -1 / 12.0L,
  -1 / 24.0L,
  -19 / 720.0L,
  -3 / 160.0L,
  -863 / 60480.0L,
  -275 / 24192.0L,
  -33953 / 3628800.0L,
  -8183 / 1036800.0L,
  -3250433 / 479001600.0L,
  -4671 / 788480.0L,
  -13695779093 / 2615348736000.0L,
  -2224234463 / 475517952000.0L,
  -132282840127 / 31384184832000.0L,
  -2639651053 / 689762304000.0L,
  -1.11956703448001e14 / 3.201186852864e16L,
  -5.0188465e7 / 1.5613165568e10L,
  -2.334028946344463e15 / 7.86014494949376e17L
};

DOUBLE *_get_moulton_coeffs(int n) {
  DOUBLE *cs = (DOUBLE *) malloc(sizeof(DOUBLE) * n);
  for (int i = 0; i < n; i++) {
    cs[i] = 0;
  }
  cs[0] = 1;
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      cs[i] += cs[j] / (DOUBLE)(i - j + 1);
    }
    cs[i] *= -1;
  }
  return cs;
}

void _backward_diff(int n, int len, DOUBLE *out) {
  out[0] = 1;
  for (int i = 1; i < len; i++) {
    if (i < n) {
      DOUBLE val = out[i - 1] * (n - i) / i;
      out[i] = powl(-1, i) * fabsl(val);
      continue;
    }
    out[i] = 0;
  }
}

void get_corrector_coeffs(int n, DOUBLE *out) {
  for (int i = 0; i < n; i++) {
    out[i] = 0;
  }
  DOUBLE *diff = (DOUBLE *) malloc(sizeof(DOUBLE) * n);
  for (int i = 0; i < n; i++) {
    _backward_diff(i + 1, n, diff);
    for (int j = 0; j < n; j++) {
      out[j] += (DOUBLE)(diff[j] * _MOULTON_COEFFS[i]);
    }
  }
  free(diff);
}

void get_predictor_coeffs(int n, DOUBLE *out) {
  for (int i = 0; i < n; i++) {
    out[i] = 0;
  }
  DOUBLE *diff = (DOUBLE *) malloc(sizeof(DOUBLE) * n);
  for (int i = 0; i < n; i++) {
    DOUBLE c = 0;
    _backward_diff(i + 1, n, diff);
    for (int j = 0; j < i + 1; j++) {
      c += _MOULTON_COEFFS[j];
    }
    for (int j = 0; j < n; j++) {
      out[j] += (DOUBLE)(diff[j] * c);
    }
  }
  free(diff);
}
