#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

long min_parallel_length =  1;

void insertionSort(int *A, long i_start, long i_end);
void merge(int *B, long i_start, long i_middle, long i_end,int *A);
void splitMerge(int *B, long i_start, long i_end, int *A);
void splitMergeP(int *B, long i_start, long i_end, int *A);


void print_v(int *v, long l) {
  printf("\n");
  for (long i = 0; i < l; i++) {
    if (i != 0 && (i % 110 == 0)) {
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

int sorted(int *v, long l){
  for(int i = 0; i< l-1; ++i){
    if(v[i] > v[i+1])
      return 0;
  }
  return 1;
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
void splitMerge(int *B, long i_start, long i_end, int *A) {
#ifdef INSERT_SORT
  if ((i_end - i_start) < 256) {
    insertionSort(A, i_start, i_end);
    return;
  }
#else
  if ((i_end - i_start) < 2)
    // A[i_start] = B[i_start];
    return;
#endif
  long i_middle = (long)((i_end + i_start) / 2);
  splitMerge(A, i_start, i_middle, B);
  splitMerge(A, i_middle, i_end, B);
  merge(B, i_start, i_middle, i_end, A);
}

// i_start is inclusive; i_end is exclusive (A[i_end] is not in the set).
void splitMergeP(int *B, long i_start, long i_end, int *A) {
  // printf("available threads: %i, %i \n",omp_get_num_threads(),omp_get_thread_num());
  long length = (i_end - i_start);

  if(length < min_parallel_length){
    splitMerge(B,i_start,i_end,A);
    return;
  }

  long i_middle = (long)((i_end + i_start) / 2);
#pragma omp task shared(A, B) firstprivate(i_start, i_middle)
  splitMergeP(A, i_start, i_middle, B);

#pragma omp task shared(A, B) firstprivate(i_end, i_middle)
  splitMergeP(A, i_middle, i_end, B);

#pragma omp taskwait
  merge(B, i_start, i_middle, i_end, A);
}

//  Left source half is A[ iBegin:iMiddle-1].
// Right source half is A[iMiddle:iEnd-1   ].
// Result is            B[ iBegin:iEnd-1   ].
void merge(int *A, long i_start, long i_middle, long i_end, int *B) {
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