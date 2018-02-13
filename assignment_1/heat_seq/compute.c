#include "compute.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

/// Replace/erase the following line:
#include "ref1.c"

#define IDX(i, j, colls) (i * colls + j)
// / printf("%E", res[IDX(i, j, M)]);
#define CROSS_RATIO 0.14644660940672627
#define DIAGONAL_RATIO 0.10355339059327377

double calcKernel(int i, int j, int rows, int colls, const Indexes *indexes) {
  double corners_sum = 0;
  double cross_sum = 0;
  int start = 0;
  int idx = IDX(i, j, colls);
  double cross = 0;
  for (int c = 1; c < 8; c += 2) {
    cross += indexes[idx].index[c];
  }

  double diag = 0;
  for (int d = 0; d < 9; d += 2) {
    diag += indexes[idx].index[d];
  }
  diag -= indexes[idx].index[4];
  // remove from sums if at each end of cylinder
  if (j == 0) {
    cross -= indexes[idx].index[3];
    diag -= indexes[idx].index[0];
    diag -= indexes[idx].index[6];
  }
  if (j == colls - 1) {
    cross -= indexes[idx].index[5];
    diag -= indexes[idx].index[2];
    diag -= indexes[idx].index[8];
  }
  return cross * CROSS_RATIO + diag * DIAGONAL_RATIO;
}

void precomputeIndexes(int cols, int rows, Indexes *indexes) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {

      int k = 0;
      for (int ii = -1; ii < 2; ++ii) {
        for (int jj = -1; jj < 2; ++jj) {
          int vert_cell_idx = i + ii;
          // wrapping
          if (vert_cell_idx < 0) {
            vert_cell_idx = rows - 1;
          }
          if (vert_cell_idx > rows - 1) {
            vert_cell_idx = 0;
          }
          indexes[IDX(i, j, cols)].index[k] = IDX(vert_cell_idx, j + jj, cols);
          ++k;
        }
      }
    }
  }
}

double minTemp(const double *res, int length) {
  double tmin = 0;
  for (int i = 0; i < length; ++i) {
    if (res[i] < tmin) {
      tmin = res[i];
    }
  }
  return tmin;
}

double maxTemp(const double *res, int length) {
  double tmax = 0;
  for (int i = 0; i < length; ++i) {
    if (res[i] > tmax) {
      tmax = res[i];
    }
  }
  return tmax;
}

double avgTemp(const double *res, int length) {
  double temps = 0;
  for (int i = 0; i < length; ++i) {
    temps += res[i];
  }
  return temps / length;
}

void do_compute(const struct parameters *p, struct results *r) {
  int M = p->M; // rows (cylinder circumference plane)
  int N = p->N; // cols (height of cylinder)
  const int length_cells = N * M;
  const int length_alloc = length_cells * sizeof(double);
  double *res = (double *)malloc(length_alloc);
  memcpy(res, p->tinit, length_alloc);

  // compute indexes
  printf("sizeof one Index struct %lli\n", sizeof(Indexes));
  Indexes *indexes = (Indexes *)malloc(length_cells * sizeof(Indexes));
  precomputeIndexes(N, M, indexes);
  // for(int i = 0;i<5;++i){
  //   for(int k = 0; k<9;++k){
  //     printf("%i  ", indexes[i].index[k]);
  //   }
  //   printf("\n%i:%i\n",i,length_cells);
  // }

  // statistics

  double tmin = minTemp(res, length_cells);
  double tmax = maxTemp(res, length_cells);
  double temp_avg = avgTemp(res, length_cells);
  double max_diff = 0;
  int iter = 0;
  for (; iter < p->maxiter; ++iter) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        int idx = IDX(i, j, M);
        double conductivity = p->conductivity[idx];
        double new_temp = calcKernel(i, j, N, M, indexes);
        double temp = conductivity * res[idx] + (1 - conductivity) * new_temp;
        double diff_temp = temp - res[idx];
        if (diff_temp > max_diff) {
          max_diff = diff_temp;
        }
        res[idx] = temp;
      }
      // printf("\n\n");
    }
    if (iter % p->period == 0) {
      r->niter = iter;
      r->tmin = tmin;
      r->tmax = tmax;
      r->tavg = temp_avg;
      r->maxdiff = max_diff;
      r->time = 10;
      report_results(p, r);
    }
    tmin = minTemp(res, length_cells);
    tmax = maxTemp(res, length_cells);
    temp_avg = avgTemp(res, length_cells);
    if (p->threshold > max_diff) {
      break;
    }
  }

  r->niter = iter;
  r->tmin = tmin;
  r->tmax = tmax;
  r->tavg = temp_avg;
  r->maxdiff = max_diff;
  r->time = 10;
  report_results(p, r);

  free(res);
  free(indexes);
}
