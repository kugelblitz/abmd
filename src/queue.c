#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "queue.h"


struct _Queue {
  double t0, h;
  int head, tail, size, dim;
  int capacity, block_size;
  DOUBLE* array;
  DOUBLE* diffs;
  DOUBLE* pdiffs;
  DOUBLE* temp;
  double *lgr_ws;
  double *lgr_ws_nolast;
  DOUBLE *lgr_nom;
};

Queue *create_queue(int capacity, int block_size) {
  int dim = block_size / 2;
  Queue* queue = (Queue *) malloc(sizeof(Queue));
  queue->t0 = 0;
  queue->h = 1;
  queue->capacity = capacity;
  queue->block_size = block_size;
  queue->dim = dim;
  queue->head = queue->size = 0;
  queue->tail = block_size * (capacity - 1);
  queue->array = (DOUBLE *) malloc(capacity * block_size * sizeof(DOUBLE));
  queue->diffs = (DOUBLE *) malloc(capacity * dim * sizeof(DOUBLE));
  queue->pdiffs = (DOUBLE *) malloc(capacity * dim * sizeof(DOUBLE));
  queue->temp = (DOUBLE *) malloc(2 * dim * sizeof(DOUBLE));
  queue->lgr_ws = (double *) malloc(capacity * sizeof(double));
  queue->lgr_ws_nolast = (double *) malloc((capacity - 1) * sizeof(double));
  queue->lgr_nom = (DOUBLE *) malloc(dim * sizeof(DOUBLE));
  return queue;
}

void destroy_queue(Queue *queue) {
  free(queue->array);
  free(queue->diffs);
  free(queue->pdiffs);
  free(queue->temp);
  free(queue->lgr_ws);
  free(queue->lgr_ws_nolast);
  free(queue->lgr_nom);
  free(queue);
}

int is_full(Queue *q) {
  return q->size == q->capacity;
}

int is_empty(Queue *q) {
  return q->size == 0;
}

int get_capacity(Queue *q) {
  return q->capacity;
}

void set_t0(Queue *q, double t0) {
  q->t0 = t0;
}

double _compute_wj(Queue *q, int j, int len) {
  /*
   * Computes j-th barycentric Lagrange weight (eq. 3.2)
  */
  double w = 1;
  double h = q->h;
  for (int k = 0; k < len; k++) {
    if (k == j) continue;
    w /= (j - k) * h;  // ts[j] - ts[k] simplified
  }
  return w;
}

void _initialize_weights(Queue *q) {
  int n = get_capacity(q);
  for (int i = 0; i < n - 1; i++) {
    q->lgr_ws[i] = _compute_wj(q, i, n);
    q->lgr_ws_nolast[i] = _compute_wj(q, i, n - 1);
  }
  q->lgr_ws[n - 1] = _compute_wj(q, n - 1, n);
}

void set_step(Queue *q, double h) {
  q->h = h;
  _initialize_weights(q);
}

int _get_t_index(Queue *q, double t, int last_known) {
  double t0 = q->t0;
  double h = q->h;
  int n = get_capacity(q);
  if (!last_known) {
    n -= 1;
  }
  double t1 = t0 + (n - 1) * h;
  int hsgn = (t1 > t0) - (t1 < t0);
  t0 *= hsgn;
  t1 *= hsgn;
  t *= hsgn;
  h *= hsgn;

  if (t0 <= t && t <= t1 && fmod((t - t0), h) < 1e-13) {
    return (int) round((t - t0) / h);
  }
  return -1;
}

void _evaluate(Queue *q, double t, int *idxs, int idxs_len, int last_known,
               int dim_start, DOUBLE *out) {

  int t_idx = _get_t_index(q, t, last_known);
  if (t_idx != -1) {
    DOUBLE *x = get(q, t_idx);
    for (int j = 0; j < idxs_len; j++) {
      int idx = (idxs == NULL) ? j : idxs[j];
      out[j] = x[idx];
    }
    return;
  }

  int n = get_capacity(q);
  double *ws = q->lgr_ws;
  if (!last_known) {
    ws = q->lgr_ws_nolast;
    n -= 1;
  }

  DOUBLE *nom = q->lgr_nom;
  memset(nom, 0, q->dim * sizeof(DOUBLE));
  DOUBLE denom = 0;
  for (int i = 0; i < n; i++) {
    double tti = t - (q->t0 + i * q->h);
    denom += ws[i] / tti;
    DOUBLE *x = &get(q, i)[dim_start];
    for (int j = 0; j < idxs_len; j++) {
      int idx = (idxs == NULL) ? j : idxs[j];
      nom[j] += ws[i] * x[idx] / tti;
    }
  }
  for (int j = 0; j < idxs_len; j++) {
    out[j] = nom[j] / denom;
  }
}

void evaluate_x_all(Queue *q, double t, DOUBLE *out) {
  _evaluate(q, t, NULL, q->dim, 1, 0, out);
}

void evaluate_x_idxs(Queue *q, double t, int *idxs, int idxs_len, DOUBLE *out) {
  _evaluate(q, t, idxs, idxs_len, 1, 0, out);
}

void evaluate_xdot(Queue *q, double t, int *idxs, int idxs_len,
                   int last_known, DOUBLE *out) {
  _evaluate(q, t, idxs, idxs_len, last_known, q->dim, out);
}

DOUBLE* get(Queue *q, int block_idx) {
  int q_size = get_capacity(q);
  if (block_idx < 0) {
    return get(q, q_size + block_idx % q_size);
  }
  int array_idx = (q->head + q->block_size * block_idx) %
                       (q->block_size * q_size);
  return &q->array[array_idx];
}

DOUBLE* push(Queue *q) {
  if (is_full(q))
    return NULL;
  q->tail = (q->tail + q->block_size) % (q->block_size * q->capacity);
  q->size = q->size + 1;
  return &(q->array[q->tail]);
}

DOUBLE* pop(Queue *q) {
  if (is_empty(q))
    return NULL;
  DOUBLE *address = &q->array[q->head];
  q->head = (q->head + q->block_size) % (q->block_size * q->capacity);
  q->size = q->size - 1;
  q->t0 += q->h;
  return address;
}

DOUBLE* peek_left(Queue* q) {
  if (is_empty(q))
    return NULL;
  return &q->array[q->head];
}

DOUBLE* peek_right(Queue* q) {
  if (is_empty(q))
    return NULL;
  return &q->array[q->tail];
}

DOUBLE* get_diff(Queue *q, int i) {
  return &q->diffs[i * q->dim];
}

DOUBLE* update_diffs(Queue *q, int predicted) {
  int qsize = q->size;
  int dim = q->dim;
  int dim_size = dim * sizeof(DOUBLE);
  DOUBLE *old = q->temp;
  DOUBLE *new = &q->temp[dim];
  memcpy(new, &peek_right(q)[dim], dim_size);

  int ndiffs = is_full(q) ? qsize : qsize + 1;

  DOUBLE *diffs = q->diffs;
  if (predicted) {
    memcpy(q->pdiffs, q->diffs, qsize * dim_size);
    diffs = q->pdiffs;
  }

  for (int i = 0; i < ndiffs - 1; i++) {
    DOUBLE *idiff = &diffs[i * dim];
    memcpy(old, idiff, dim_size);
    memcpy(idiff, new, dim_size);
    for (int j = 0; j < dim; j++) {
      new[j] -= old[j];
    }
  }

  memcpy(&diffs[(ndiffs - 1) * dim], new, dim_size);
  return diffs;
}


int test() {
  Queue* queue = create_queue(4, 2);
  DOUBLE *address;
  address = push(queue);
  address[0] = 10;
  address[1] = 15;

  assert(peek_left(queue) == peek_right(queue));
  assert(peek_left(queue)[0] == 10);
  assert(peek_left(queue)[1] == 15);

  address = push(queue);
  address[0] = 20;
  address[1] = 25;

  assert(peek_left(queue) != peek_right(queue));
  assert(peek_left(queue)[0] == 10);
  assert(peek_right(queue)[0] == 20);
  assert(peek_right(queue)[1] == 25);

  address = push(queue);
  address[0] = 30;
  address[1] = 35;
  address = push(queue);
  address[0] = 40;
  address[1] = 45;

  assert(pop(queue)[0] == 10);
  assert(peek_left(queue)[1] == 25);
  assert(peek_right(queue)[0] == 40);

  assert(get(queue, 0)[0] == 20);
  assert(get(queue, 1)[1] == 35);
  assert(get(queue, 2)[0] == 40);
  assert(get(queue, -3)[0] == 30);
  assert(get(queue, -4)[0] == 20);
  assert(get(queue, -7)[1] == 35);

  printf("Queue tests passed");
  return 0;
}

int main_() {
  int qsize = 3;
  int dim = 2;
  Queue *q = create_queue(qsize, dim * 2);
  set_step(q, -1);
  set_t0(q, 0);
  DOUBLE *one = push(q);
  DOUBLE *two = push(q);
  DOUBLE *three = push(q);
  one[0] = 0;
  one[1] = 16;
  two[0] = 1;
  two[1] = 9;
  three[0] = 4;
  three[1] = 4;

  DOUBLE *out = (DOUBLE *) malloc(dim * sizeof(DOUBLE));
  evaluate_x_all(q, -1.5, out);
  printf("%Le, %Le\n", out[0], out[1]);

  pop(q);
  DOUBLE *next = push(q);
  next[0] = 9;
  next[1] = 1;
  evaluate_x_all(q, -4, out);
  printf("%Le, %Le\n", out[0], out[1]);
  
  return 0;
}
