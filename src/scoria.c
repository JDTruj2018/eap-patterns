#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

void print_tsum(double val, int32_t ndim) {
  printf("TSum: %.2f, ndim: %" PRId32 "\n", val, ndim);
}

void print_face_num(int32_t * fnum, const int len) {
  for (int i = 0; i < len; ++i) 
    printf("i: %d\t%" PRId32 "\n", i, fnum[i]);
  printf("\n");
}

void scoria_read_1(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1) {
  for (size_t i = 0; i < N1; ++i) {
    res[i] = buffer[ind1[i] - 1];
  }
}

void scoria_write_1(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1) {
  for (size_t i = 0; i < N1; ++i) {
    buffer[ind1[i] - 1] = input[i];
  }
}