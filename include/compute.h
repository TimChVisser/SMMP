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

#define VISUALIZATION

#define CROSS_RATIO 0.14644660940672627
#define CROSS_RATIO_MISS 0.1098349609834961
#define DIAGONAL_RATIO 0.10355339059327377
#define DIAGONAL_RATIO_MISS 0.0517766952965



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
