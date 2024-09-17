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

void scoria_read_1_many_2(double *resA, double *resB, const double *bufferA, const double *bufferB, const size_t N, const size_t *ind1, const size_t N1) {
  for (size_t i = 0; i < N1; ++i) {
    resA[i] = bufferA[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resB[i] = bufferB[ind1[i] - 1];
  }
}

void scoria_read_1_many_3(double *resA, double *resB, double* resC, const double *bufferA, const double *bufferB, const double *bufferC, const size_t N, const size_t *ind1, const size_t N1) {
  for (size_t i = 0; i < N1; ++i) {
    resA[i] = bufferA[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resB[i] = bufferB[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resC[i] = bufferC[ind1[i] - 1];
  }
}

void scoria_read_1_many_4(double *resA, double *resB, double* resC, double* resD, const double *bufferA, const double *bufferB, const double *bufferC, const double *bufferD, const size_t N, const size_t *ind1, const size_t N1) {
  for (size_t i = 0; i < N1; ++i) {
    resA[i] = bufferA[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resB[i] = bufferB[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resC[i] = bufferC[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resD[i] = bufferD[ind1[i] - 1];
  }
}

void scoria_read_1_many_6(double *resA, double *resB, double* resC, double* resD, double* resE, double* resF, const double *bufferA, const double *bufferB, const double *bufferC, double *bufferD, const double *bufferE, const double *bufferF, const size_t N, const size_t *ind1, const size_t N1) {
  for (size_t i = 0; i < N1; ++i) {
    resA[i] = bufferA[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resB[i] = bufferB[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resC[i] = bufferC[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resD[i] = bufferD[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resE[i] = bufferE[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resF[i] = bufferF[ind1[i] - 1];
  }
}

void scoria_read_1_many_7(double *resA, double *resB, double* resC, double* resD, double* resE, double* resF, double* resG, const double *bufferA, const double *bufferB, const double *bufferC, const double *bufferD, const double *bufferE, const double *bufferF, const double *bufferG, const size_t N, const size_t *ind1, const size_t N1) {
  for (size_t i = 0; i < N1; ++i) {
    resA[i] = bufferA[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resB[i] = bufferB[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resC[i] = bufferC[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resD[i] = bufferD[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resE[i] = bufferE[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resF[i] = bufferF[ind1[i] - 1];
  }
  for (size_t i = 0; i < N1; ++i) {
    resG[i] = bufferG[ind1[i] - 1];
  }
}
