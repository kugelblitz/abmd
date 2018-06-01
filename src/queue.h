#ifndef DDE_QUEUE_H
#define DDE_QUEUE_H

#include "abm.h"

typedef struct _Queue Queue;

Queue *create_queue(int capacity, int block_size);
void destroy_queue(Queue *queue);
int is_full(Queue *q);
int is_empty(Queue *q);
int get_capacity(Queue *q);
DOUBLE* get(Queue *q, int block_idx);
DOUBLE* push(Queue *q);
DOUBLE* pop(Queue *q);
DOUBLE* peek_left(Queue* q);
DOUBLE* peek_right(Queue* q);
DOUBLE* get_diffs_r(Queue *q);
DOUBLE* get_diffs_w(Queue *q);
DOUBLE* get_last_diff(Queue *q);
void swap_diffs(Queue *q);
void update_diffs(Queue* q);
void set_t0(Queue *q, double t0);
void set_step(Queue *q, double h);
void evaluate_x_all(Queue *q, double t, DOUBLE *out);
void evaluate_x_idxs(Queue *q, double t, int *idxs, int idxs_len, DOUBLE *out);
void evaluate_xdot(Queue *q, double t, int *idxs, int idxs_len,
                   int last_known, DOUBLE *out);

#endif //DDE_QUEUE_H
