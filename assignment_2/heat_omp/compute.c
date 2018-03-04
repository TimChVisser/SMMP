#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#include <omp.h>
#include <stdio.h>

/// Replace/erase the following line:
#include "ref1.c"

void do_compute(const struct parameters* p, struct results *r)
{
//omp_set_num_threads(p->nthreads);
#include "ref2.c"
}
