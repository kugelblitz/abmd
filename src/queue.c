#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "queue.h"


struct _Queue {
  double t0, h;
  int head, tail, size, dim;
  int capacity, block_size;
  DOUBLE* array;
  DOUBLE* diffs;
  DOUBLE* pdiffs;
  DOUBLE* temp;
  double *_lgr_ws, *_lgr_denom;
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
  queue->_lgr_ws = (double *) malloc(capacity * sizeof(double));
  queue->_lgr_denom = (double *) malloc(block_size * sizeof(double));
  return queue;
}

void destroy_queue(Queue *queue) {
  free(queue->array);
  free(queue->diffs);
  free(queue->pdiffs);
  free(queue->temp);
  free(queue->_lgr_ws);
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

double get_t0(Queue *q) {
  return q->t0;
}

void set_step_size(Queue *q, double h) {
  q->h = h;
}

double get_step_size(Queue *q) {
  return q->h;
}

double _compute_wj(Queue *q, int j) {
  /*
   * Computes j-th barycentric Lagrange weight (eq. 3.2)
  */
  double w = 1;
  int n = get_capacity(q);
  double h = get_step_size(q);
  for (int k = 0; k < n; k++) {
    if (k == j) continue;
    w /= (j - k) * h;  // ts[j] - ts[k] simplified
  }
  return w;
}

void _initialize_ws(Queue *q) {
  int n = get_capacity(q);
  for (int i = 0; i < n; i++) {
    q->_lgr_ws[i] = _compute_wj(q, i);
  }
}

void evaluate(Queue *q, double t, double *out) {
  int n = get_capacity(q);
  int dim = q->block_size;
  double t0 = get_t0(q);
  double h = get_step_size(q);
  DOUBLE *nom = malloc(dim * sizeof(DOUBLE));
  DOUBLE denom = 0;
  for (int j = 0; j < n; j++) {
    denom += q->_lgr_ws[j] / (t - t0 - j * h);
    DOUBLE *block = get(q, j);
    for (int i = 0; i < dim; i++) {
      nom[i] += q->_lgr_ws[j] * block[i] / (t - t0 - j * h);
    }
  }
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
