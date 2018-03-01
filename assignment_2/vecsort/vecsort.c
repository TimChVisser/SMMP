
#include "../mergesort/sort.h"
#include <ctype.h>
#include <getopt.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Ordering of the vector */
typedef enum Ordering { ASCENDING, DESCENDING, RANDOM } Order;

#define MAX_INER_LENGTH 100000

int debug = 0;
int seed = 42;
int cores1 = 2;
int cores2 = 3;
int iterations = 100;

struct Array {
  long size;
  int *data;
};

long random_at_most(long max) {
  unsigned long
      // max <= RAND_MAX < ULONG_MAX, so this is okay.
      num_bins = (unsigned long)max + 1,
      num_rand = (unsigned long)RAND_MAX + 1, bin_size = num_rand / num_bins,
      defect = num_rand % num_bins;

  long x;
  do {
    x = rand();
  }
  // This is carefully written not to overflow
  while (num_rand - defect <= (unsigned long)x);

  // Truncated division is intentional
  return x / bin_size;
}

void createInputVectors(struct Array *A, long length, long max_iner_length,
                        Order order) {
  srand(seed);
  for (long i = 0; i < length; ++i) {
    A[i].size = random_at_most(max_iner_length - 1) + 1;
    A[i].data = (int *)malloc(A[i].size * sizeof(int));

    if (A[i].data == NULL) {
      fprintf(stderr, "Malloc failed...\n");
      return;
    }
  }
  // populate with numbers
  for (long i = 0; i < length; ++i) {
    long max = A[i].size;
    switch (order) {
    case ASCENDING:
      for (long j = 0; j < A[i].size; j++) {
        A[i].data[j] = (int)j;
      }
      break;
    case DESCENDING:
      for (long j = 0; j < A[i].size; ++j) {
        A[i].data[j] = (int)(max - j);
      }
      break;
    case RANDOM:
      for (long j = 0; j < A[i].size; j++) {
        A[i].data[j] = rand();
      }
      break;
    }
  }
}

void vecsort(long length, long max_iner_length, Order order) {

  double time_iner = 0;
  double time_total = 0;

  struct Array *main_arr =
      (struct Array *)malloc(length * sizeof(struct Array));

  if (main_arr == NULL) {
    fprintf(stderr, "Malloc failed - main_arr...\n");
    return;
  }

  int **scratch_arr = (int **)malloc(length * sizeof(int *));
  for (long i = 0; i < cores1; ++i) {
    scratch_arr[i] = (int *)malloc(max_iner_length * sizeof(int));
  }

  for (long i = 0; i < iterations; ++i) {

    createInputVectors(main_arr, length, max_iner_length, order);

    int thread_id = 0;
    struct timeval tv3, tv4;
    double time = 0;
    gettimeofday(&tv3, NULL);

#pragma omp parallel for num_threads(cores1) schedule(guided)                  \
    shared(scratch_arr, main_arr) reduction(+ : time)
    struct timeval tv1, tv2;
    for (int i = 0; i < length; ++i) {

      memcpy(scratch_arr[thread_id], main_arr[i].data,
             main_arr[i].size * sizeof(int));

      gettimeofday(&tv1, NULL);

#pragma omp parallel num_threads(cores2)
      {
#pragma omp single nowait
        splitMergeP(scratch_arr[thread_id], 0, main_arr[i].size,
                    main_arr[i].data);
      }
      gettimeofday(&tv2, NULL);

      time = time + ((double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
                     (double)(tv2.tv_sec - tv1.tv_sec));
    }

    gettimeofday(&tv4, NULL);
    time_total += (double)(tv4.tv_usec - tv3.tv_usec) / 1000000 +
                  (double)(tv4.tv_sec - tv3.tv_sec);
    time_iner += time;

    // print result
    if (debug == 1) {
      for (long i = 0; i < length; ++i) {
        print_v(main_arr[i].data, main_arr[i].size);
      }
    }

    // deallocate iner vectors
    for (long i = 0; i < length; ++i) {
      free(main_arr[i].data);
    }
  }
  // average times
  printf("%i,%i,%i,%i,%e\n", cores1, cores2, length, max_iner_length,
         time_total / iterations);

  free(main_arr);
  for (long i = 0; i < cores1; ++i) {
    free(scratch_arr[i]);
  }
  free(scratch_arr);
}

int main(int argc, char **argv) {

  int c;

  long length = 1e4;
  long max_iner_length = 1e5;

  Order order = DESCENDING;

  /* Read command-line options. */
  while ((c = getopt(argc, argv, "adrgl:s:v:c:m:i:L:")) != -1) {
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
    case 'm':
      max_iner_length = atol(optarg);
      break;
    case 'g':
      debug = 1;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'c':
      cores1 = atoi(optarg);
      break;
    case 'v':
      cores2 = atoi(optarg);
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

  // if(debug) {
  //   print_v(/* ... */);
  // }

  // /* Sort */
  vecsort(length, max_iner_length, order);

  // if(debug) {
  //   print_v(/* ... */);
  // }

  return 0;
}
