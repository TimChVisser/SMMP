//#define OPENMP
//#define PTHREAD

size_t i, j;

/* alias input parameters */
const double (*restrict tinit)[p->N]
                              [p->M] = (const double (*)[p->N][p->M])p->tinit;
const double (*restrict cinit)[p->N][p->M] = (const double (*)[p->N][p->M])
                                                 p->conductivity;

/* allocate grid data */
const size_t h = p->N + 2;
const size_t w = p->M + 2;
double (*restrict g1)[h][w] = malloc(h * w * sizeof(double));
double (*restrict g2)[h][w] = malloc(h * w * sizeof(double));

/* allocate halo for conductivities */
double (*restrict c)[h][w] = malloc(h * w * sizeof(double));

struct timeval before;


/* set initial temperatures and conductivities */
//#pragma omp parallel for collapse(2)
for(i = 1; i < h - 1; ++i){
  for (j = 1; j < w - 1; ++j) {
    (*g1)[i][j] = (*tinit)[i - 1][j - 1];
    (*c)[i][j] = (*cinit)[i - 1][j - 1];
  }
}
/* smear outermost row to border */
//#pragma omp parallel for
for (j = 1; j < w - 1; ++j) {
  (*g1)[0][j] = (*g2)[0][j] = (*g1)[1][j];
  (*g1)[h - 1][j] = (*g2)[h - 1][j] = (*g1)[h - 2][j];
}

/* compute */
size_t iter;
double (*restrict src)[h][w] = g2;
double (*restrict dst)[h][w] = g1;

/*
 * If initialization should be included in the timings
 * could be a point of discussion.
 */
gettimeofday(&before, NULL);
for (iter = 1; iter <= p->maxiter; ++iter) {
#ifdef GEN_PICTURES
  do_draw(p, iter, h, w, src);
#endif
  /* swap source and destination */
  {
    void *tmp = src;
    src = dst;
    dst = tmp;
  }

  /* initialize halo on source */
  do_copy(h, w, src);

  double maxdiff = 0.0;
/* compute */

#ifdef PTHREAD
   printf("hello1%i\n",omp_get_num_threads() );
  pthread_t thread_ids[threads];
  struct Params params[threads];
  double results[threads];
  pthread_attr_t attributes;
  pthread_attr_init(&attributes);
  pthread_attr_setdetachstate(&attributes ,  PTHREAD_CREATE_JOINABLE );
   pthread_attr_setscope( &attributes ,  PTHREAD_SCOPE_SYSTEM );
  for(size_t i =0; i < threads; ++i){
    params[i].src = src;
    params[i].dst = dst;
    params[i].c = c;
    params[i].start = i * ((h-1) / threads);
    if(threads-1 == i)
    {
        params[i].end = h - 1;
    }else{
        params[i].end = (i + 1) * ((h-1) / threads) -1;
    }
    params[i].maxdiff = &(results[i]);
    pthread_create(&thread_ids[i],&attributes,&performeCalc,&params[i]);
  }
  for(size_t i =0 ;i < threads; ++i){
    pthread_join(thread_ids[i],NULL);
  }
  pthread_attr_destroy(&attributes);
  // reduction of maxdiff
  maxdiff = results[0];
  for(size_t i =1; i < threads; ++i){
    if ((results[i]) > maxdiff){
        maxdiff = (results[i]);
    }
  }

#else
#ifdef OPENMP
#pragma omp parallel num_threads(threads)
{
#pragma omp for
  for (i = 1; i < h - 1; ++i) {
    for (j = 1; j < w - 1; ++j) {

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
}
#else
  for (i = 1; i < h - 1; ++i) {
    for (j = 1; j < w - 1; ++j) {

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
#endif // OpenMP
#endif // PTHREAD
  r->maxdiff = maxdiff;
  if (maxdiff < p->threshold) {
    iter++;
    break;
  }
  /* conditional reporting */
  if (iter % p->period == 0) {
    fill_report(p, r, h, w, dst, src, iter, &before);
    if (p->printreports)
      report_results(p, r);
  }
}

/* report at end in all cases */
iter--;
fill_report(p, r, h, w, dst, src, iter, &before);

free(c);
free(g2);
free(g1);
