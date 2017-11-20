#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "type.h"
#include "coeffs.h"
#include "poly.h"
#include "plot.h"


#define ABM_ORDER 11


typedef struct {
    void (*call)(DOUBLE, DOUBLE[], void*, DOUBLE*);
    int dim;
    int ndelays;
    DOUBLE *delays;
    void *context;
} FContext;

void copy(DOUBLE* from, DOUBLE* to, int count) {
    for (int i = 0; i < count; i++) {
        to[i] = from[i];
    }
}

void predict(int i, DOUBLE h, DOUBLE *sol, DOUBLE *rhss, int dim, DOUBLE *out,
             DOUBLE *coeffs) {
    for (int j = 0; j < dim; j++) {
        out[j] = sol[(i - 1) * dim + j];
    }
    for (int j = 0; j < ABM_ORDER; j++) {
        DOUBLE *rhs = &rhss[(i - j - 1) * dim];
        DOUBLE c = coeffs[j];
        for (int k = 0; k < dim; k++) {
            out[k] += c * h * rhs[k];
        }
    }
}

void correct(int i, DOUBLE h, DOUBLE *sol, DOUBLE *rhss, int dim, DOUBLE *out,
             DOUBLE *coeffs) {
    for (int j = 0; j < dim; j++) {
        out[j] = sol[(i - 1) * dim + j];
    }

    for (int j = 0; j < ABM_ORDER; j++) {
        DOUBLE *rhs = &rhss[(i - j) * dim];
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

void rk4_step(void (*f)(DOUBLE, DOUBLE[], void*, DOUBLE*), DOUBLE h, DOUBLE t,
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
        copy(state, &in1[i * dim], dim);
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

void compute_func(DOUBLE t, DOUBLE *state, void *fcontext, DOUBLE *out) {

    FContext *fctx = (FContext *)fcontext;
    int dim = fctx->dim;
    DOUBLE *states = malloc(sizeof(DOUBLE) * fctx->ndelays * dim);
    for (int i = 0; i < fctx->ndelays; i++) {
        DOUBLE delay = fctx->delays[i];
        if (delay == 0) {
            copy(state, &states[i * dim], dim);
            continue;
        }
        rk4_step(fctx->call, -delay, t, state, dim, fctx->ndelays, fctx->context, &states[i * dim]);
    }
    fctx->call(t, states, fctx->context, out);
    free(states);
}


DOUBLE *get_delayed_states(int i, DOUBLE *ts, DOUBLE *sol, DOUBLE h, int interpolation_order,
                           int extrapolation_order, DOUBLE delays[], int dim, int delays_number) {
    DOUBLE *states = malloc(sizeof(DOUBLE) * dim * delays_number);
    for (int j = 0; j < delays_number; j++) {
        DOUBLE delay = delays[j];
        DOUBLE t_delayed = ts[i] - delay;
        if (delay >= 0) {
            if (!fmod(delay / h, 1)) {
                int row = i - (int)(delay / h);
                for (int ii = 0; ii < dim; ii++) {
                    states[j * dim + ii] = sol[row * dim + ii];
                }
            } else {
                // Interpolation
                int points_number = interpolation_order + 1;
                int right_i = i - (int)(delay / h);
                if ((int)(delay / h) > 0) {
                    right_i += 1;
                }
                DOUBLE *xs = &ts[right_i - points_number + 1];
                DOUBLE *ys = &sol[(right_i - points_number + 1) * dim];
                lagrange(t_delayed, xs, ys, dim, points_number, &states[j * dim]);
            }
        } else {
            // Extrapolation
            int points_number = extrapolation_order + 1;
            DOUBLE *xs = &ts[i - points_number + 1];
            DOUBLE *ys = &sol[(i - points_number + 1) * dim];
            lagrange(t_delayed, xs, ys, dim, points_number, &states[j * dim]);
        }
    }
    return states;
}


DOUBLE *abm(void (*f)(DOUBLE t, DOUBLE states[], void *context, DOUBLE *out), DOUBLE t0, DOUBLE t1, DOUBLE h,
            DOUBLE init[], DOUBLE delays[], int dim, int delays_number, void *context) {
/*  Integrates the DDE system given by `f` using Adams—Bashforth—Moulton method of order `ABM_ORDER`.
 *
 *  Integration is performed from `t0` to `t1` with constant step `h`,
 *  `init` being an array of length `dim` containing the system's state at `t0`.
 *
 *  f arguments:
 *      t - point in time (variable which is being integrated with respect to);
 *      states - array of length `dim * delays_number` storing system states corresponding to `delays`;
 *      context - argument being passed to `f` through the `abm` function;
 *      out - function output (`dim` number of DOUBLEs).
 */
    FContext fcontext;
    fcontext.call = f;
    fcontext.context = context;
    fcontext.dim = dim;
    fcontext.ndelays = delays_number;
    fcontext.delays = delays;

    DOUBLE predictor_coeffs[ABM_ORDER];
    DOUBLE corrector_coeffs[ABM_ORDER];
    get_predictor_coeffs(ABM_ORDER, predictor_coeffs);
    get_corrector_coeffs(ABM_ORDER, corrector_coeffs);

    int RK_STEPS_IN_ABM = 4;

    int interpolation_order = 4;
    int extrapolation_order = 1;

    DOUBLE max_positive_delay = 0;
    DOUBLE max_negative_delay = 0;
    for (int i = 0; i < delays_number; i++) {
        if (delays[i] > max_positive_delay) {
            max_positive_delay = delays[i];
        }
        if (delays[i] < max_negative_delay) {
            max_negative_delay = delays[i];
        }
    }

    int n = (int)(1 + (t1 - t0) / h);
    int extra_steps = fmod(max_positive_delay, h) ?
                      (int)(max_positive_delay / h) + 1 : (int)(max_positive_delay / h);

    extra_steps += interpolation_order - 1;
//    printf("Max positive delay: %Lf\n", max_positive_delay);
//    printf("Extra steps: %d\n", extra_steps);
    DOUBLE rk4_t1 = t0 + (ABM_ORDER - 1 + extra_steps) * h;
//    printf("rk4_t1: %Lf\n", rk4_t1);
    DOUBLE rk4_h = h / (DOUBLE)RK_STEPS_IN_ABM;
    int rk4_n = (int)(1 + (rk4_t1 - t0) / rk4_h);

    DOUBLE *ts = (DOUBLE *)malloc(sizeof(DOUBLE) * n);
    DOUBLE *rhss = (DOUBLE *)malloc(sizeof(DOUBLE) * n * dim);
    DOUBLE *sol = (DOUBLE *)malloc(sizeof(DOUBLE) * n * dim);
    DOUBLE *rk4_sol = (DOUBLE *)malloc(sizeof(DOUBLE) * rk4_n * dim);

    // Setting initial conditions for RK4 solution
    copy(init, rk4_sol, dim);

    for (int i = 0; i < n; i++) {
        ts[i] = t0 + i * h;
    }


    // Doing rk4_n RK4 steps
    for (int i = 1; i < rk4_n; i++) {
        DOUBLE t = t0 + rk4_h * i;
        rk4_step(compute_func, rk4_h, t, &rk4_sol[(i - 1) * dim], dim, 1, &fcontext, &rk4_sol[i * dim]);
    }

    // Writing data from RK4 to the solution array
    int k = 0;
    for (int i = 0; i < rk4_n; i += RK_STEPS_IN_ABM) {
        copy(&rk4_sol[i * dim], &sol[k * dim], dim);
        k++;
    }
    free(rk4_sol);

    // Initializing right-hand sides
    for (int i = extra_steps; i < k; i++) {
        DOUBLE *states = get_delayed_states(i, ts, sol, h, interpolation_order, extrapolation_order,
                                            delays, dim, delays_number);
        f(ts[i], states, context, &rhss[i * dim]);
        free(states);
    }


    // If ABM_ORDER = 2: rk4_t1 = 2000, t = 4000
    int start_index = (int)(rk4_t1 / h) + 1;

    for (int i = start_index; i < n; i++) {

        predict(i, h, sol, rhss, dim, &sol[i * dim], predictor_coeffs);
        DOUBLE *states = get_delayed_states(i, ts, sol, h, interpolation_order,
                                            extrapolation_order, delays, dim, delays_number);
        f(ts[i], states, context, &rhss[i * dim]);
        free(states);

        correct(i, h, sol, rhss, dim, &sol[i * dim], corrector_coeffs);
        states = get_delayed_states(i, ts, sol, h, interpolation_order,
                                    extrapolation_order, delays, dim, delays_number);
        f(ts[i], states, context, &rhss[i * dim]);
        free(states);
    }
    free(ts);
    free(rhss);
    return sol;
}

#ifdef DEBUG

void calc_difference(void f(DOUBLE, DOUBLE*, void*, DOUBLE*)) {
    DOUBLE init[] = {-3844e5, 0, 0, 1023};
    DOUBLE t0 = 0;
    DOUBLE t1 = 3e7;
    DOUBLE h = 2000;
    DOUBLE delay = 3000;
    int dim = 4;
    int n = (int)(1 + (t1 - t0) / h);
    DOUBLE delays[2] = {0, delay};
    DOUBLE *sol = abm(f, t0, t1, h, init, delays, dim, 2, NULL);
    DOUBLE *sol_reversed = malloc(sizeof(DOUBLE) * n * dim);

    for (int i = 0; i < n; i++) {
        copy(&sol[(n - i - 1) * dim], &sol_reversed[i * dim], dim);
        sol_reversed[i * dim + 1] *= -1;
        sol_reversed[i * dim + 3] *= -1;
    }

    DOUBLE delays_[2] = {0, -delay};
    DOUBLE init_back[] = {sol[(n - 1) * dim],
                          -sol[(n - 1) * dim + 1],
                          sol[(n - 1) * dim + 2],
                          -sol[(n - 1) * dim + 3]};
    DOUBLE *sol_back = abm(f, t0, t1, h, init_back, delays_, dim, 2, NULL);

    DOUBLE *diff = malloc(sizeof(DOUBLE) * dim * n);

    for (int i = 0; i < n; i++) {
        DOUBLE x1 = sol_reversed[i * dim];
        DOUBLE x2 = sol_back[i * dim];
        DOUBLE y1 = sol_reversed[i * dim + 2];
        DOUBLE y2 = sol_back[i * dim + 2];
        diff[i * dim] = t0 + i * h;
        diff[i * dim + 2] = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
    }

    plot("diff.png", diff, dim, n);
    free(sol);
    free(sol_reversed);
    free(sol_back);
    free(diff);
}

#endif
