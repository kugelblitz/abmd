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
  DOUBLE* _array;
  DOUBLE **xarray, **dxarray;
  DOUBLE *x_backup;
  DOUBLE* diffs_r;
  DOUBLE* diffs_w;
  DOUBLE* last_diff;
  double *lgr_ws;
  double *lgr_ws_nolast;
  DOUBLE *lgr_nom;
};

Queue *create_queue(int capacity, int dim) {
  int block_size = 2 * dim;
  Queue* queue = (Queue *) malloc(sizeof(Queue));
  queue->t0 = 0;
  queue->h = 1;
  queue->capacity = capacity;
  queue->block_size = block_size;
  queue->dim = dim;
  queue->head = queue->size = 0;
  queue->tail = capacity - 1;
  queue->_array = (DOUBLE *) malloc(capacity * block_size * sizeof(DOUBLE));
  queue->xarray = (DOUBLE **) malloc(capacity * sizeof(DOUBLE *));
  queue->dxarray = (DOUBLE **) malloc(capacity * sizeof(DOUBLE *));
  queue->x_backup = (DOUBLE *) malloc(dim * sizeof(DOUBLE));

  for (int i = 0; i < capacity; i++) {
    queue->xarray[i] = &queue->_array[i * block_size];
    queue->dxarray[i] = &queue->_array[i * block_size + dim];
  }

  queue->diffs_r = (DOUBLE *) malloc((capacity - 1) * dim * sizeof(DOUBLE));
  queue->diffs_w = (DOUBLE *) malloc((capacity - 1) * dim * sizeof(DOUBLE));
  queue->last_diff = (DOUBLE *) malloc(dim * sizeof(DOUBLE));
  queue->lgr_ws = (double *) malloc(capacity * sizeof(double));
  queue->lgr_ws_nolast = (double *) malloc((capacity - 1) * sizeof(double));
  queue->lgr_nom = (DOUBLE *) malloc(dim * sizeof(DOUBLE));
  return queue;
}

void destroy_queue(Queue *queue) {
  free(queue->_array);
  free(queue->xarray);
  free(queue->dxarray);
  free(queue->x_backup);
  free(queue->diffs_r);
  free(queue->diffs_w);
  free(queue->last_diff);
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

DOUBLE* get_x(Queue *q, int block_idx) {
  return q->xarray[(q->head + block_idx) % q->capacity];
}

DOUBLE* get_dx(Queue *q, int block_idx) {
  return q->dxarray[(q->head + block_idx) % q->capacity];
}

void backup_last_x(Queue *q) {
  q->xarray[q->tail] = q->x_backup;
}

void restore_last_x(Queue *q) {
  int i = (q->tail * q->block_size) % (q->capacity * q->block_size);
  q->xarray[q->tail] = &q->_array[i];
}

void _evaluate(Queue *q, double t, int *idxs, int idxs_len, int last_known,
               DOUBLE *(*get)(Queue *, int), DOUBLE *out) {

  int t_idx = _get_t_index(q, t, last_known);
  if (t_idx != -1) {
    DOUBLE *x = get(q, t_idx);
    if (idxs == NULL)
      memcpy(out, x, idxs_len * sizeof(DOUBLE));
    else
      for (int j = 0; j < idxs_len; j++) {
        out[j] = x[idxs[j]];
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
    double coef = ws[i] / tti;
    denom += coef;
    DOUBLE *x = get(q, i);
    if (idxs == NULL) {
      for (int j = 0; j < idxs_len; j++) {
        nom[j] += coef * x[j];
      }
    } else {
      for (int j = 0; j < idxs_len; j++) {
        nom[j] += coef * x[idxs[j]];
      }
    }
  }
  for (int j = 0; j < idxs_len; j++) {
    out[j] = nom[j] / denom;
  }
}

void evaluate_x_all(Queue *q, double t, DOUBLE *out) {
  _evaluate(q, t, NULL, q->dim, 1, get_x, out);
}

void evaluate_x_idxs(Queue *q, double t, int *idxs, int idxs_len, DOUBLE *out) {
  _evaluate(q, t, idxs, idxs_len, 1, get_x, out);
}

void evaluate_xdot(Queue *q, double t, int *idxs, int idxs_len,
                   int last_known, DOUBLE *out) {
  _evaluate(q, t, idxs, idxs_len, last_known, get_dx, out);
}

DOUBLE* push(Queue *q) {
  if (is_full(q))
    return NULL;
  q->tail = (q->tail + 1) % q->capacity;
  q->size += 1;
  return q->xarray[q->tail];
}

DOUBLE* pop(Queue *q) {
  if (is_empty(q))
    return NULL;
  DOUBLE *address = q->xarray[q->head];
  q->head = (q->head + 1) % q->capacity;
  q->size -= 1;
  q->t0 += q->h;
  return address;
}

DOUBLE* peek_left(Queue* q) {
  if (is_empty(q))
    return NULL;
  return q->xarray[q->head];
}

DOUBLE* peek_right_x(Queue* q) {
  if (is_empty(q))
    return NULL;
  return q->xarray[q->tail];
}

DOUBLE* peek_right_dx(Queue* q) {
  if (is_empty(q))
    return NULL;
  return q->dxarray[q->tail];
}

DOUBLE* get_diffs_r(Queue *q) {
  return q->diffs_r;
}

DOUBLE* get_diffs_w(Queue *q) {
  return q->diffs_w;
}

DOUBLE* get_last_diff(Queue *q) {
  return q->last_diff;
}

void swap_diffs(Queue *q) {
  DOUBLE *diffs_r = q->diffs_r;
  q->diffs_r = q->diffs_w;
  q->diffs_w = diffs_r;
}

void update_diffs(Queue *q) {
  int size = q->size;
  int diffs_len = q->capacity - 1;
  int dim = q->dim;

  DOUBLE *right = peek_right_dx(q);
  DOUBLE *diffs_r = q->diffs_r;
  DOUBLE *diffs_w = q->diffs_w;

  if (!is_full(q)) {
    for (int i = 0; i < dim; i++) {
      DOUBLE new = *right++;
      for (int j = 0; j < size - 1; j++) {
        *diffs_w++ = new;
        new -= *diffs_r++;
      }
      *diffs_w = new;
      diffs_w += diffs_len - size + 1;
      diffs_r += diffs_len - size + 1;
    }
    return;
  }

  DOUBLE *last_diff = q->last_diff;

  for (int i = 0; i < dim; i++) {
    DOUBLE new = *right++;
    for (int j = 0; j < size - 1; j++) {
      *diffs_w++ = new;
      new -= *diffs_r++;
    }
    *last_diff++ = new;
  }
}


int test() {
  Queue* queue = create_queue(4, 2);
  DOUBLE *address;
  address = push(queue);
  address[0] = 10;
  address[1] = 15;

  assert(peek_left(queue) == peek_right_x(queue));
  assert(peek_left(queue)[0] == 10);
  assert(peek_left(queue)[1] == 15);

  address = push(queue);
  address[0] = 20;
  address[1] = 25;

  assert(peek_left(queue) != peek_right_x(queue));
  assert(peek_left(queue)[0] == 10);
  assert(peek_right_x(queue)[0] == 20);
  assert(peek_right_x(queue)[1] == 25);

  address = push(queue);
  address[0] = 30;
  address[1] = 35;
  address = push(queue);
  address[0] = 40;
  address[1] = 45;

  assert(pop(queue)[0] == 10);
  assert(peek_left(queue)[1] == 25);
  assert(peek_right_x(queue)[0] == 40);

  assert(get_x(queue, 0)[0] == 20);
  assert(get_x(queue, 1)[1] == 35);
  assert(get_x(queue, 2)[0] == 40);

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
