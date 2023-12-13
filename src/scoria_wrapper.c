#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#include "kernels.h"

void print_tsum(double val, int32_t ndim) {
  printf("TSum: %.2f, ndim: %" PRId32 "\n", val, ndim);
}

void print_face_num(int32_t * fnum, const int len) {
  for (int i = 0; i < len; ++i) 
    printf("i: %d\t%" PRId32 "\n", i, fnum[i]);
  printf("\n");
}

i_type get_intrin(const int intrin) {
  i_type intrinsic;
  if (intrin == 0) intrinsic = NONE;
  else if (intrin == 1) intrinsic == AVX;
  else intrinsic == SVE;

  return intrinsic;
}

void read_single_thread_0_wrapper(double *res, const double *buffer, const size_t N, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  read_single_thread_0(res, buffer, N, intrinsic);
}

void read_single_thread_1_wrapper(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  read_single_thread_1(res, buffer, N, ind1, intrinsic);
}

void read_single_thread_2_wrapper(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1, const size_t *ind2, const size_t N2, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  read_single_thread_2(res, buffer, N, ind1, ind2, intrinsic);
}

void read_multi_thread_0_wrapper(double *res, const double *buffer, const size_t N, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  read_multi_thread_0(res, buffer, N, nthreads, intrinsic);
}

void read_multi_thread_1_wrapper(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  read_multi_thread_1(res, buffer, N, ind1, nthreads, intrinsic);
}

void read_multi_thread_2_wrapper(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1, const size_t *ind2, const size_t N2, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  read_multi_thread_2(res, buffer, N, ind1, ind2, nthreads, intrinsic);
}


void write_single_thread_0_wrapper(double *buffer, const double *input, const size_t N, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  write_single_thread_0(buffer, input, N, intrinsic);
}

void write_single_thread_1_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  write_single_thread_1(buffer, input, N, ind1, intrinsic);
}

void write_single_thread_2_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const size_t *ind2, const size_t N2, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  write_single_thread_2(buffer, input, N, ind1, ind2, intrinsic);
}

void write_multi_thread_0_wrapper(double *buffer, const double *input, const size_t N, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  write_multi_thread_0(buffer, input, N, nthreads, intrinsic);
}

void write_multi_thread_1_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  write_multi_thread_1(buffer, input, N, ind1, nthreads, intrinsic);
}

void write_multi_thread_2_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const size_t *ind2, const size_t N2, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  write_multi_thread_2(buffer, input, N, ind1, ind2, nthreads, intrinsic);
}


void writeadd_single_thread_0_wrapper(double *buffer, const double *input, const size_t N, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  writeadd_single_thread_0(buffer, input, N, intrinsic);
}

void writeadd_single_thread_1_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  writeadd_single_thread_1(buffer, input, N, ind1, intrinsic);
}

void writeadd_single_thread_2_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const size_t *ind2, const size_t N2, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  writeadd_single_thread_2(buffer, input, N, ind1, ind2, intrinsic);
}

void writeadd_multi_thread_0_wrapper(double *buffer, const double *input, const size_t N, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  writeadd_multi_thread_0(buffer, input, N, nthreads, intrinsic);
}

void writeadd_multi_thread_1_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  writeadd_multi_thread_1(buffer, input, N, ind1, nthreads, intrinsic);
}

void writeadd_multi_thread_2_wrapper(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const size_t *ind2, const size_t N2, const size_t nthreads, const int intrin) {
  i_type intrinsic = get_intrin(intrin);

  writeadd_multi_thread_2(buffer, input, N, ind1, ind2, nthreads, intrinsic);
}
