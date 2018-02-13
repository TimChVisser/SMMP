#include "compute.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

/// Replace/erase the following line:
#include "ref1.c"

#define IDX(i, j, cols) (i * cols + j)
// / printf("%E", res[IDX(i, j, M)]);
#define CROSS_RATIO 0.14644660940672627
#define DIAGONAL_RATIO 0.10355339059327377
// #define VISUALIZATION

double calcKernel(int i, int j, int rows, int colls, const Indexes *indexes,
                  const double *temps) {
  double corners_sum = 0;
  double cross_sum = 0;
  int start = 0;
  int idx = IDX(i, j, colls);
  double cross = 0;
  for (int c = 1; c < 8; c += 2) {
    cross += temps[indexes[idx].index[c]];
  }

  double diag = 0;
  for (int d = 0; d < 9; d += 2) {
    diag += temps[indexes[idx].index[d]];
  }
  diag -= temps[indexes[idx].index[4]];
  // remove from sums if at each end of cylinder
  if (j == 0) {
    // printf("left %i,%i:%i \n",i,j,indexes[IDX(i,j,colls)].index[3]);
    cross -= temps[indexes[idx].index[3]];
    diag -= temps[indexes[idx].index[0]];
    diag -= temps[indexes[idx].index[6]];
  }
  if (j == colls - 1) {
    // printf("right %i,%i:%i \n",i,j,indexes[IDX(i,j,colls)].index[5]);
    cross -= temps[indexes[idx].index[5]];
    diag -= temps[indexes[idx].index[2]];
    diag -= temps[indexes[idx].index[8]];
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

void sendStatistics(const struct parameters *p, struct results *stats,
                    const double *temps, const double *last_temps, int length,
                    int iteration, double time, double *max_diff) {

  double tmin = minTemp(temps, length);
  double tmax = maxTemp(temps, length);
  double temp_avg = avgTemp(temps, length);

  for (size_t i = 0; i < length; ++i) {
    double diff = temps[i] - last_temps[i];
    if (diff > *max_diff) {
      (*max_diff) = diff;
    }
  }
  stats->niter = iteration;
  stats->tmin = tmin;
  stats->tmax = tmax;
  stats->tavg = temp_avg;
  stats->maxdiff = (*max_diff);
  stats->time = time;
  report_results(p, stats);
}

void do_compute(const struct parameters *p, struct results *r) {
  int cols = p->M; // (height of cylinder)
  int rows = p->N; //(cylinder circumference plane)
  const int length_cells = cols * rows;
  const int length_alloc = length_cells * sizeof(double);
  double *res = (double *)malloc(length_alloc);
  double *last_res = (double *)malloc(length_alloc);
  memcpy(res, p->tinit, length_alloc);
  memcpy(last_res, res, length_alloc);

  // compute indexes
  printf("sizeof one Index struct %lli\n", sizeof(Indexes));
  Indexes *indexes = (Indexes *)malloc(length_cells * sizeof(Indexes));
  precomputeIndexes(cols, rows, indexes);
  // for(int i = 0;i<5;++i){
  //   for(int k = 0; k<9;++k){
  //     printf("%i  ", indexes[i].index[k]);
  //   }
  //   printf("\n%i:%i\n",i,length_cells);
  // }

  // statistics
  // for (size_t i = 0; i < 9; ++i) {
  //   if (i % 3 == 0) {
  //     printf("\n");
  //   }
  //   printf("%i  ", indexes[IDX(1, 101, cols)].index[i]);
  // }
  double start = 0; //(double)time(NULL);
  int iter = 0;
  double time = 0;
  for (; iter < p->maxiter; ++iter) {

#ifdef VISUALIZATION
    begin_picture(iter, cols, rows, p->io_tmin, p->io_tmax);
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        draw_point(i, j, res[IDX(i, j, cols)]);
      }
    }
    end_picture();
#endif

    memcpy(last_res, res, length_alloc);

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        int idx = IDX(i, j, cols);
        double conductivity = p->conductivity[idx];
        // printf("%e\n",conductivity );
        double new_temp = calcKernel(i, j, rows, cols, indexes, last_res);
        double temp =
            conductivity * last_res[idx] + (1 - conductivity) * new_temp;
        // double temp = new_temp;
        res[idx] = temp;
      }
      // printf("\n\n");
    }
    // double end = (double)time(NULL);
    time = 0;//end - start;
    if (iter % p->period == 0) {
      double max_diff = 0;
      sendStatistics(p, r, res, last_res, length_cells, iter, time, &max_diff);
      if (p->threshold > max_diff) {
        break;
      }
    }
  }

  double max_diff = 0;
  sendStatistics(p, r, res, last_res, length_cells, iter, time, &max_diff);
  free(res);
  free(last_res);
  free(indexes);
}
