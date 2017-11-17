#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "type.h"

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
    DOUBLE *cs = malloc(sizeof(DOUBLE) * n);
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
    DOUBLE *diff = malloc(sizeof(DOUBLE) * n);
    for (int i = 0; i < n; i++) {
        _backward_diff(i + 1, n, diff);
        for (int j = 0; j < n; j++) {
            out[j] += (double)(diff[j] * _MOULTON_COEFFS[i]);
        }
    }
    free(diff);
}

void get_predictor_coeffs(int n, DOUBLE *out) {
    for (int i = 0; i < n; i++) {
        out[i] = 0;
    }
    DOUBLE *diff = malloc(sizeof(DOUBLE) * n);
    for (int i = 0; i < n; i++) {
        DOUBLE c = 0;
        _backward_diff(i + 1, n, diff);
        for (int j = 0; j < i + 1; j++) {
            c += _MOULTON_COEFFS[j];
        }
        for (int j = 0; j < n; j++) {
            out[j] += (double)(diff[j] * c);
        }
    }
    free(diff);
}