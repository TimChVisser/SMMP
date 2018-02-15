#include "compute.h"
#include "immintrin.h"

///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// //
// __m256d _mm256_set_pd (double e3, double e2, double e1, double e0) // //
// __m256d _mm256_shuffle_pd (__m256d a, __m256d b, const int imm8)   // //
// __m256d _mm256_hadd_pd (__m256d a, __m256d b)                      // //
// __m256d _mm256_permute2f128_pd (__m256d a, __m256d b, int imm8)
// __m256d _mm256_permute_pd (__m256d a, int imm8)
// __m256d _mm256_blend_pd (__m256d a, __m256d b, const int imm8)

//////////////////////////////////////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////

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

double maxDiff(const double *temps, const double *last_temps, double length) {
  double max_diff = 0;
  for (size_t i = 0; i < length; ++i) {
    double diff = abs(temps[i] - last_temps[i]);
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
  stats->maxdiff = maxDiff(temps, last_temps, length);
  stats->time = time;
  report_results(p, stats);
}

void visualize(const struct parameters *p, const double *res, int iter,
               int rows, int cols) {
  begin_picture(iter, cols, rows, p->io_tmin, p->io_tmax);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 1; j < cols-1; ++j) {
      draw_point(i, j, res[IDX(i, j, cols)]);
    }
  }
  end_picture();
}

// int cache_size = sizeof(double) * 4 int row_size =
//     cahce_size * 10; // ?? select based on table enlargement so every row
//     starts
//                      // on the beginning of cache line

// __m256d singleVal(__m256d val1, __m256d val2) {
//   __m256d val1_prep = prepSingleVal(val1);
//   __m256d val2_prep = prepSingleVal(val2);
//   return _mm256_permute2f128_pd(val1_prep, val2_prep,
//                                 0b011000) // x_4 x_3 x_2 x_1
// }

// __m256d prepSingleVal(__m256d val) {
//   __m256d swaped_left = _mm256_permute2f128_pd(val, val, 0b011011);
//   return _mm256_shuffle_pd(val, swaped_left, 0b1011); // R R x_2 x_1
// }

// __m256d prepTwoVal(__m256d val) {
//   __m256d swaped_left = _mm256_permute2f128_pd(val, val, 0b011011);
//   return _mm256_shuffle_pd(val, swaped_left, 0b0100); // R R x_2 x_1
// }

// __m256d twoVals(__m256d val1, __m256d val2, __m256d val3) {
//   __m256d x1ox3 = _mm256_permute2f128_pd(val1, val2, 0b010000);
//   __m256d x1x3x3x4 = _mm256_shuffle_pd(x1ox3, val2, 0);
//   __m256d x2ox2o = _mm256_permute2f128_pd(val1, val1, 0b011011);
//   __m256d x2x3x3x4 = _mm256_blend_pd(x1x3x3x4, x2ox2o, 0b0001);

//   __m256d x1x2oo = prepTwoVal(val1);
//   __m256d x4x5oo = preptwoVal(val3);
//   __m256d x1x2x4x5 = _mm256_permute2f128_pd(x1x2oo, x4x5oo, 0b011000)

//       __m256d sum = _mm256_hadd_pd(x1x2x4x5, x2x3x3x4);
//   return __m256d _mm256_permute_pd(sum, 0b0110); //// x_1 x_2 x_3 x_4
// }

// // this is used for calculation of single elements in cross on top and bottom
// __m256d top1 _mm256_load_pd(last_temps);
// __m256d top2 _mm256_load_pd(last_temps + cache_size);
// __m256d bottom1 _mm256_load_pd(last_temps + 2 * row_size);
// __m256d bottom2 _mm256_load_pd(last_temps + 2 * row_size + cahce_size);

// // bottom and top cross of 4 fields
// __m256d top_cross = singleVal(top1, top2);
// __m256d bottom_cross = singleVal(bottom1, bottom2);
// __m256d cross = _mm256_add_pd(top_cross, bottom_cross);

// // middle part of 4 crosses
// __m256d fields_reg[3];
// for (int i = 0; i < 3; ++i) {
//   fields_reg[i] = _mm256_load_pd(last_temps + 1 * row_size + i * cache_size);
// }
// // sum of all elements in cross
// cross =
//     _mm256_add_pd(cross, twoVals(fields_reg[0], fields_reg[1],
//     fields_reg[2]));

// // diagonal
// __m256d diagonal = _mm256_setzero_pd();
// for (int parts = 0; parts < 3; parts += 2) {
//   for (int i = 0; i < 3; ++i) {
//     fields_reg[i] =
//         _mm256_load_pd(last_temps + parts * row_size + i * cache_size);
//   }
//   // sum of all elements in cross
//   diagonal = _mm256_add_pd(
//       diagonal, twoVals(fields_reg[0], fields_reg[1], fields_reg[2]));
// }

__m256d twoVals(__m256d val1) {
  __m256d swaped = _mm256_permute2f128_pd(val1, val1, 0b00001);
  return _mm256_add_pd(val1, swaped);
}

__m256d prepValToMiddle(__m256d val1) {
  __m256d swaped = _mm256_permute2f128_pd(val1, val1, 0b1001);
  return _mm256_shuffle_pd(val1, swaped, 0b0001);
}

__m256d calcKernels(__m256d top, __m256d middle, __m256d bottom,
                    __m256d conductivity, double diag_ratio,
                    double cross_ratio) {
  conductivity = prepValToMiddle(conductivity);

  __m256d diagonal = _mm256_add_pd(twoVals(top), twoVals(bottom));
  __m256d cross = _mm256_add_pd(twoVals(middle), prepValToMiddle(top));
  cross = _mm256_add_pd(cross, prepValToMiddle(top));
  __m256d diagonal_weighted = _mm256_mul_pd(cross, _mm256_set1_pd(diag_ratio));
  __m256d cross_weighted = _mm256_mul_pd(cross, _mm256_set1_pd(cross_ratio));
  __m256d sum = _mm256_add_pd(diagonal_weighted, cross_weighted);

  __m256d cond_minus = _mm256_sub_pd(_mm256_set1_pd(1), conductivity);
  __m256d result = _mm256_mul_pd(cond_minus, sum);
  // __m256d old_res = _mm256_mul_pd(conductivity, prepValToMiddle(middle));
  return result;//_mm256_add_pd(old_res, result);
}
void avxDataManipulCompute(const struct parameters *p, double *temps,
                           const double *last_temps, int rows, int cols,
                           int r_cols) {
  // skip first and last row -> these rows will be done separately
  for (int i = 1; i < rows - 1; ++i) {
    // first column
    // __m256d top = _mm256_loadu_pd(last_temps + IDX(i - 1, 0, r_cols));
    // __m256d middle = _mm256_loadu_pd(last_temps + IDX(i, 0, r_cols));
    // __m256d bottom = _mm256_loadu_pd(last_temps + IDX(i + 1, 0, r_cols));
    // __m256d conductivity = _mm256_loadu_pd(p->conductivity + IDX(i, 0, cols));
    // __m256d total = calcKernels(top, middle, bottom, conductivity,
    //                             DIAGONAL_RATIO_MISS, CROSS_RATIO_MISS);

    // _mm256_storeu_pd(temps + IDX(i, 0, r_cols), total);

    for (int j = 3; j < cols - 1; j += 2) {
      __m256d top = _mm256_loadu_pd(last_temps + IDX(i - 1, j - 1, r_cols));
      __m256d middle = _mm256_loadu_pd(last_temps + IDX(i, j - 1, r_cols));
      __m256d bottom = _mm256_loadu_pd(last_temps + IDX(i + 1, j - 1, r_cols));
      __m256d conductivity =
          _mm256_loadu_pd(p->conductivity + IDX(i, j - 1, cols));
      __m256d total = calcKernels(top, middle, bottom, conductivity,
                                DIAGONAL_RATIO, CROSS_RATIO);
      _mm256_storeu_pd(temps + IDX(i, j, r_cols), total);
    }
    // last column
    // top = _mm256_loadu_pd(last_temps + IDX(i - 1, cols-2, r_cols));
    // middle = _mm256_loadu_pd(last_temps + IDX(i, 0, cols-2));
    // bottom = _mm256_loadu_pd(last_temps + IDX(i + 1, cols-2, r_cols));
    // conductivity = _mm256_loadu_pd(p->conductivity + IDX(i,cols-2, cols));
    // total = calcKernels(top, middle, bottom, conductivity,
    //                             DIAGONAL_RATIO_MISS, CROSS_RATIO_MISS);

    // _mm256_storeu_pd(temps + IDX(i, cols-2, r_cols), total);
  }
}

void do_compute(const struct parameters *p, struct results *r) {
  int cols = p->M; // (height of cylinder)
  int rows = p->N; //(cylinder circumference plane)
  int r_cols = cols + 2;
  int r_rows = rows;

  const int length_cells = cols * rows;
  const int length_alloc = r_cols * r_rows * sizeof(double);
  const int length_alloc_s = length_cells * sizeof(double);
  double *res = (double *)malloc(length_alloc);
  double *last_res = (double *)malloc(length_alloc);
  // first and last column has zeros. p->tinit is in between
  for (size_t i = 0; i < rows;++i) {
    res[IDX(i, 0, r_cols)] = 0;
    for (size_t j = 1; j < cols;++j) {
      res[IDX(i, j, r_cols)] = p->tinit[IDX(i, j, cols)];
    }
    res[IDX(i, cols, r_cols)] = 0;
  }
  // memcpy(res+1, p->tinit, length_alloc_s);
  memcpy(last_res, res, length_alloc);

  // run main loop
  int iter = 0;
  double time = 0;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (; iter < p->maxiter; ++iter) {

#ifdef VISUALIZATION
    visualize(p, res, iter, rows, r_cols);
#endif
    memcpy(last_res, res, length_alloc);
    avxDataManipulCompute(p, res, last_res, rows, cols, r_cols);

    if (p->threshold > maxDiff(res, last_res, length_cells)) {
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

  // visualize(p, res, 0, rows, r_cols);

  gettimeofday(&tv2, NULL);
  time = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double)(tv2.tv_sec - tv1.tv_sec);
  sendStatistics(p, r, res, last_res, length_cells, iter, time);

  free(res);
  free(last_res);
}
