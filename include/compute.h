#ifndef COMPUTE_H
#define COMPUTE_H

#include "input.h"
#include "output.h"

typedef struct {
  int index[9];
} Indexes;
void precomputeIndexes(int cols, int rows, Indexes *indexes);
void do_compute(const struct parameters *p, struct results *r);
#endif
