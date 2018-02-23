
#include <ctype.h>
#include <getopt.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#define INSERT_SORT
/* Ordering of the vector */
typedef enum Ordering { ASCENDING, DESCENDING, RANDOM } Order;

int debug = 0;
int multicore = 0;
int cores = 0;
int iterations = 1;

void insertionSort(int *A, long i_start, long i_end);
void merge(int *B, long i_start, long i_middle, long i_end, const int *A);
void splitMerge(int *B, long i_start, long i_end, const int *A);
void splitMergeP(int *B, long i_start, long i_end, const int *A);
void msort(int *data, long length, Order order);

void print_v(int *v, long l) {
  printf("\n");
  for (long i = 0; i < l; i++) {
    if (i != 0 && (i % 10 == 0)) {
      printf("\n");
    }
    printf("%d ", v[i]);
  }
  printf("\n");
}

void printRange(const int *v, long start, long end) {
  printf("\n");
  for (long i = start; i < end; i++) {
    if (i != 0 && (i % 10 == 0)) {
      printf("\n");
    }
    printf("%d ", v[i]);
  }
  printf("\n");
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
int main(int argc, char **argv) {

  int c;
  int seed = 42;
  long length = 1e4;
  Order order = ASCENDING;
  int *vector;

  /* Read command-line options. */
  while ((c = getopt(argc, argv, "adrgpl:s:c:i:")) != -1) {
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

  if (debug) {
    print_v(vector, length);
  }

  free(vector);
  return 0;
}
////////////////////////////////////////////////////////

void insertionSort(int *A, long i_start, long i_end) {
  long i, j;
  int val;
  // int * start_adr = A + i_start;
  int *current_addr_j = A + i_start;
  for (i = i_start + 1; i < i_end; ++i) {
    val = *(A + i);
    j = i - 1;
    current_addr_j = A + i - 1;
    while (j >= i_start && *(current_addr_j) > val) {

      *(current_addr_j + 1) = *(current_addr_j);
      *(current_addr_j) = val;
      --current_addr_j;
      --j;
    }
  }
}

///////////

// i_start is inclusive; i_end is exclusive (A[i_end] is not in the set).
void splitMerge(int *B, long i_start, long i_end, const int *A) {
#ifdef INSERT_SORT
  if ((i_end - i_start) < 1000) {
    insertionSort(B, i_start, i_end);
    return;
  }
#else
  if ((i_end - i_start) < 2)
    return;
#endif
  long i_middle = (long)((i_end + i_start) / 2);

  splitMerge(B, i_start, i_middle, A);
  splitMerge(B, i_middle, i_end, A);
  merge(B, i_start, i_middle, i_end, A);
}

// i_start is inclusive; i_end is exclusive (A[i_end] is not in the set).
void splitMergeP(int *B, long i_start, long i_end, const int *A) {
#ifdef INSERT_SORT
  if ((i_end - i_start) < 1000) {
    insertionSort(B, i_start, i_end);
    return;
  }
#else
  if ((i_end - i_start) < 2)
    return;
#endif
  long i_middle = (long)((i_end + i_start) / 2);
#pragma omp tasks shared(A, B) firstprivate(i_start, i_middle)
  splitMergeP(B, i_start, i_middle, A);

#pragma omp tasks shared(A, B) firstprivate(i_end, i_middle)
  splitMergeP(B, i_middle, i_end, A);

#pragma omp taskwait
  merge(B, i_start, i_middle, i_end, A);
}

//  Left source half is A[ iBegin:iMiddle-1].
// Right source half is A[iMiddle:iEnd-1   ].
// Result is            B[ iBegin:iEnd-1   ].
void merge(int *B, long i_start, long i_middle, long i_end, const int *A) {
  long i = i_start;
  long j = i_middle;
  for (long k = i_start; k < i_end; k++) {
    // If left run head exists and is <= existing right run head.
    if (i < i_middle && (j >= i_end || A[i] <= A[j])) {
      B[k] = A[i];
      i = i + 1;
    } else {
      B[k] = A[j];
      j = j + 1;
    }
  }
}
//////////////////////// PARALEL

////////////////////////////////////////////
void msort(int *data, long length, Order order) {
  int *B = (int *)malloc(length * sizeof(int));
  if (B == NULL) {
    fprintf(stderr, "Malloc failed - B...\n");
    return;
  }
  int *A = data;
  // print_v(A+10, length-10);
  // insertionSort(A, 10, length);
  // printf("---------------------\n");
  // print_v(A+10, length-10);
  double time;

  for (long i = 0; i < iterations; ++i) {
    initVector(order, A, length);

    memcpy(B, A, length);
    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);
    if (multicore) {
#pragma omp parallel shared(B, A), firstprivate(length)
      {
#pragma omp single
        splitMerge(B, 0, length, A);
      }
    } else {
      splitMerge(B, 0, length, A);
    }
    gettimeofday(&tv2, NULL);
    time += (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
            (double)(tv2.tv_sec - tv1.tv_sec);
  }
  // average times
  printf("%e\n", time / iterations);

  data = B;
  free(B);
}
