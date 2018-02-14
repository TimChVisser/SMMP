#ifndef COMPUTE_H
#define COMPUTE_H

#include "input.h"
#include "output.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#define IDX(i, j, cols) (i * cols + j)
// / printf("%E", res[IDX(i, j, M)]);
#define CROSS_RATIO 0.14644660940672627
#define DIAGONAL_RATIO 0.10355339059327377
// #define VISUALIZATION



typedef struct {
  int index[9]; // indexes in temperature array. index[0] - top left cell of kernel
  				// index[8] - bottom right cell of kernel
  				// left and right border just overflows
} Indexes;

double calcKernel(int i, int j, int rows, int colls, const Indexes *indexes,
                  const double *temps);

void precomputeIndexes(int cols, int rows, Indexes *indexes);

double minTemp(const double *res, int length);
double maxTemp(const double *res, int length);
double avgTemp(const double *res, int length);
double maxDiff(const double *temps,const double *last_temps,double length);
void sendStatistics(const struct parameters *p, struct results *stats,
                    const double *temps, const double *last_temps, int length,
                    int iteration, double time);
void visualize(const struct parameters *p, const double *res, int iter, int rows,
               int cols);


void do_compute(const struct parameters *p, struct results *r);

#endif
