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
DOUBLE* get_diff(Queue *q, int i);
DOUBLE* update_diffs(Queue* q, int predicted);

#endif //DDE_QUEUE_H
