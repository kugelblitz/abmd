#ifndef DDE_QUEUE_H
#define DDE_QUEUE_H

#include "abmd.h"

typedef struct _Queue Queue;

Queue *create_queue(int capacity, int block_size);
void destroy_queue(Queue *queue);
int is_full(Queue *q);
int is_empty(Queue *q);
int get_capacity(Queue *q);
void backup_last_x(Queue *q);
void restore_last_x(Queue *q);
ABMD_DOUBLE* push(Queue *q);
ABMD_DOUBLE* pop(Queue *q);
ABMD_DOUBLE* peek_left(Queue* q);
ABMD_DOUBLE* peek_right_x(Queue* q);
ABMD_DOUBLE* peek_right_dx(Queue* q);
ABMD_DOUBLE* get_diffs_r(Queue *q);
ABMD_DOUBLE* get_diffs_w(Queue *q);
ABMD_DOUBLE* get_last_diff(Queue *q);
void swap_diffs(Queue *q);
void update_diffs(Queue* q);
void set_t0(Queue *q, double t0);
void set_step(Queue *q, double h);
void qset_delays_poly_degree(Queue *q, int deg);
void qset_pointsave_poly_degree(Queue *q, int deg);
void evaluate_x_all(Queue *q, double t, ABMD_DOUBLE *out);
void evaluate_x_idxs(Queue *q, double t, int *idxs, int idxs_len, ABMD_DOUBLE *out);
void evaluate_dx(Queue *q, double t, int *idxs, int idxs_len,
                 int last_known, ABMD_DOUBLE *out);

#endif //DDE_QUEUE_H
