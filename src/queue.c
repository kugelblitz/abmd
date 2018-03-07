#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "queue.h"


struct _Queue {
  int head, tail, size;
  int capacity, block_size;
  DOUBLE* array;
};

Queue *create_queue(int capacity, int block_size) {
  Queue* queue = (Queue *)malloc(sizeof(Queue));
  queue->capacity = capacity;
  queue->block_size = block_size;
  queue->head = queue->size = 0;
  queue->tail = block_size * (capacity - 1);
  queue->array = (DOUBLE *)malloc(capacity * block_size * sizeof(DOUBLE));
  return queue;
}

void destroy_queue(Queue *queue) {
  free(queue->array);
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
