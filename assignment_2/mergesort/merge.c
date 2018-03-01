
#include <ctype.h>
#include <getopt.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

// #define INSERT_SORT
#include "sort.h"

/* Ordering of the vector */
typedef enum Ordering { ASCENDING, DESCENDING, RANDOM } Order;
void msort(int *data, long length, Order order);
void initVector(Order order, int *vector, long length);

int debug = 0;
int multicore = 0;
int cores = 1;
int iterations = 1;


int main(int argc, char **argv) {

  int c;
  int seed = 42;
  long length = 1e4;
  Order order = ASCENDING;
  int *vector;

  /* Read command-line options. */
  while ((c = getopt(argc, argv, "adrgpl:s:c:i:L:")) != -1) {
    switch (c) {
    case 'a':
      order = ASCENDING;
      break;
    case 'd':
      order = DESCENDING;
      break;
    case 'r':
      order = RANDOM;
      break;
    case 'l':
      length = atol(optarg);
      break;
    case 'g':
      debug = 1;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'p':
      multicore = 1;
      break;
    case 'c':
      cores = atoi(optarg);
      break;
    case 'i':
      iterations = atoi(optarg);
      break;
    case 'L':
      min_parallel_length = atol(optarg);
      break;
    case '?':
      if (optopt == 'l' || optopt == 's') {
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      } else if (isprint(optopt)) {
        fprintf(stderr, "Unknown option '-%c'.\n", optopt);
      } else {
        fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
      }
      return -1;
    default:
      return -1;
    }
  }

  /* Seed such that we can always reproduce the same random vector */
  srand(seed);

  /* Allocate vector. */
  vector = (int *)malloc(length * sizeof(int));
  if (vector == NULL) {
    fprintf(stderr, "Malloc failed...\n");
    return -1;
  }

  /* Fill vector. */
  initVector(order, vector, length);

  if (debug) {
    print_v(vector, length);
  }

  /* Sort */
  msort(vector, length, order);

  // if (debug) {
  //   print_v(vector, length);
  // }

  // free(vector);
  return 0;
}

void initVector(Order order, int *vector, long length) {
  switch (order) {
  case ASCENDING:
    for (long i = 0; i < length; i++) {
      vector[i] = (int)i;
    }
    break;
  case DESCENDING:
    for (long i = 0; i < length; i++) {
      vector[i] = (int)(length - i);
    }
    break;
  case RANDOM:
    for (long i = 0; i < length; i++) {
      vector[i] = rand();
    }
    break;
  }
}

void msort(int *data, long length, Order order) {
  // int *B = (int *)malloc(length * sizeof(int));
  // if (B == NULL) {
  //   fprintf(stderr, "Malloc failed - B...\n");
  //   return;
  // }
  // int *A = data;
  // print_v(A+10, length-10);
  // insertionSort(A, 10, length);
  // printf("---------------------\n");
  // print_v(A+10, length-10);
  double time;

  for (long i = 0; i < iterations; ++i) {
    int *B = (int *)malloc(length * sizeof(int));
    if (B == NULL) {
      fprintf(stderr, "Malloc failed - B...\n");
      return;
    }

    int *A = (int *)malloc(length * sizeof(int));
    if (A == NULL) {
      fprintf(stderr, "Malloc failed - A...\n");
      return;
    }
    initVector(order, A, length);
    memcpy(B, A, length * sizeof(int));

    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);
    if (multicore) {
#pragma omp parallel num_threads(cores), shared(B, A), firstprivate(length),firstprivate(min_parallel_length)
      {
#pragma omp single
        splitMergeP(B, 0, length, A);
      }
    } else {
      splitMerge(B, 0, length, A);
    }
    gettimeofday(&tv2, NULL);
    time += (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
            (double)(tv2.tv_sec - tv1.tv_sec);

    memcpy(data,A,length * sizeof(int));
    free(A);
    free(B);
  }

  // print_v(data, length);
  // average times
  //printf("%i, %e\n", sorted(data,length),time / iterations);
  if(!multicore)
    cores = 0;
   printf("%i %i %i %i %d %e\n", min_parallel_length,multicore, cores, length, order, time / iterations);
  // data = B;
  // free(B);
}
