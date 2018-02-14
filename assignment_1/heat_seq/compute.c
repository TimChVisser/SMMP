#include "compute.h"




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

double maxDiff(const double *temps,const double *last_temps,double length){
  double max_diff = 0;
  for (size_t i = 0; i < length; ++i) {
    double diff = temps[i] - last_temps[i];
    if (diff > max_diff) {
      max_diff = diff;
    }
  }
  return max_diff;
}

void sendStatistics(const struct parameters *p, struct results *stats,
                    const double *temps, const double *last_temps, int length,
                    int iteration, double time) {

  double tmin = minTemp(temps, length);
  double tmax = maxTemp(temps, length);
  double temp_avg = avgTemp(temps, length);
  stats->niter = iteration;
  stats->tmin = tmin;
  stats->tmax = tmax;
  stats->tavg = temp_avg;
  stats->maxdiff = maxDiff(temps,last_temps,length);
  stats->time = time;
  report_results(p, stats);
}

void visualize(const struct parameters *p, const double *res, int iter, int rows,
               int cols) {
  begin_picture(iter, cols, rows, p->io_tmin, p->io_tmax);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      draw_point(i, j, res[IDX(i, j, cols)]);
    }
  }
  end_picture();
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
  printf("sizeof one Index struct %lui\n", sizeof(Indexes));
  Indexes *indexes = (Indexes *)malloc(length_cells * sizeof(Indexes));
  precomputeIndexes(cols, rows, indexes);

  // run main loop
  int iter = 0;
  double time = 0;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (; iter < p->maxiter; ++iter) {

#ifdef VISUALIZATION
    visualize(p,res,iter,rows,cols);
#endif
    memcpy(last_res, res, length_alloc);

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        int idx = IDX(i, j, cols);
        double conductivity = p->conductivity[idx];

        double new_temp = calcKernel(i, j, rows, cols, indexes, last_res);
        double temp =
            conductivity * last_res[idx] + (1 - conductivity) * new_temp;

        res[idx] = temp;
      }
    }

    if (p->threshold > maxDiff(res,last_res,length_cells)) {
        break;
    }

    if (iter % p->period == 0) {
      gettimeofday(&tv2, NULL);
      time = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
              (double)(tv2.tv_sec - tv1.tv_sec);
      sendStatistics(p, r, res, last_res, length_cells, iter, time);

      gettimeofday(&tv1, NULL);
    }
  }

  visualize(p,res,0,rows,cols);

  gettimeofday(&tv2, NULL);
      time = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
              (double)(tv2.tv_sec - tv1.tv_sec);
  sendStatistics(p, r, res, last_res, length_cells, iter, time);

  free(res);
  free(last_res);
  free(indexes);
}
