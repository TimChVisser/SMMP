#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#include <omp.h>
#include <stdio.h>
#include <pthread.h>

#include "ref1.c"

#define M_SQRT2    1.41421356237309504880
static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

struct Params{
    void *src;
    void * dst;
    void * c;
    double * maxdiff;
    size_t start;
    size_t end;
    size_t h;
    size_t w;
};

void * performeCalc(void * param) {
  struct Params * p = (struct Params *) param;
  const size_t w  = p->w;
  const size_t h = p->h;
  double (*restrict src)[h][w] =  p->src;
  double (*restrict dst)[h][w] =  p->dst;
  double (*restrict c)[h][w] =  p->c;
  size_t start = p->start;
  size_t end =p->end;

  double maxdiff = 0;
  for (size_t i = start; i < end; ++i) {
    for (size_t j = 1; j < w - 1; ++j) {

      double v = (*c)[i][j];
      double restw = 1.0 - v;

      (*dst)[i][j] = v * (*src)[i][j] +

                     ((*src)[i + 1][j] + (*src)[i - 1][j] + (*src)[i][j + 1] +
                      (*src)[i][j - 1]) *
                         (restw * c_cdir) +

                     ((*src)[i - 1][j - 1] + (*src)[i - 1][j + 1] +
                      (*src)[i + 1][j - 1] + (*src)[i + 1][j + 1]) *
                         (restw * c_cdiag);

      double diff = fabs((*dst)[i][j] - (*src)[i][j]);
      if (diff > maxdiff)
        maxdiff = diff;
    }
  }
  *(p->maxdiff) = maxdiff;
  return NULL;
}

void do_compute(const struct parameters* p, struct results *r)
{
	size_t threads = p->nthreads;
	#include "ref2.c"

 }
