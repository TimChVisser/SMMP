#include "compute.h"



#define FPOPS_PER_POINT_PER_ITERATION (                 \
        1     /* current point 1 mul */ +               \
        3 + 1 /* direct neighbors 3 adds + 1 mul */ +   \
        3 + 1 /* diagonal neighbors 3 adds + 1 mul */ + \
        2     /* final add */ +                         \
        1     /* difference old/new */                  \
        )

int main(int argc, char **argv)
{
    struct parameters p;
    struct results r;

    read_parameters(&p, argc, argv);

    do_compute(&p, &r);

    report_results(&p, &r);

    printf("%zu %zu %zu %zu %.6e %.6e\n",
           r.niter,
           p.N,
           p.M,
           p.nthreads,
           r.time,
           (double)p.N * (double)p.M *
           (double)(r.niter * FPOPS_PER_POINT_PER_ITERATION +
                    (double)r.niter / p.period) / r.time);

    return 0;
}
