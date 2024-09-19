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

// Option 1 - would not work
// Fortran -> 1-indexed
// Ranged Read as shown in Prodigy paper (except the Prodigy paper uses a summation since it is summing values of leaves)
// https://tnm.engin.umich.edu/wp-content/uploads/sites/353/2021/01/2021.02.Prodigy_HPCA2021_Camera_Ready.pdf
// ind1_orig = [1, 2, 3, 6, 7]
// buffer = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
// res_orig = [1.0, 2.0, 3.0, 6.0, 7.0]

// ind1_transformed = [1, 3, 6, 7]
// buffer = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
// res_transformed = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] - How can we skip 4.0 and 5.0? -> length can't be known before hand (in case of repeated indices in ind1, otherwise length = max(ind1_transformed))
void scoria_read_ranged_1_a(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1) {
  size_t count = 0;
  for (size_t i = 0; i < N1; ++i) {
    for(size_t j = ind1[i] - 1; j < ind1[i+1] - 1; ++j, ++count) {
      res[count] = buffer[j]; // In prodigy, this is res[i] += buffer[j] and the res array has a known length of N1
    }
  }
}

// Option 2 - could work by identifying areas of stride 1 in the indirection array ind1
// Fortran -> 1-indexed
// Ranged Read where the ind1 array contains start of offset and R array contains length of streaming access
// ind1_orig = [1, 2, 3, 6, 7]
// buffer = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
// res_orig = [1.0, 2.0, 3.0, 6.0, 7.0]

// ind1_transformed = [1, 6]
// R = [3, 2] -> length of result is 3+2 = 5
// buffer = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
// res_transformed = [1.0, 2.0, 3.0, 6.0, 7.0] -> length can be calculated by summing the values of the R array
void scoria_read_ranged_1_b(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1, const size_t *R) {
  size_t count = 0;
  for (size_t i = 0; i < N1; ++i) {
    for(size_t j = ind1[i] - 1; j < ind1[i] + R[i] - 1; ++j, ++count) {
      res[count] = buffer[j];
    }
  }
}

// Option 3 - could work by merging arrays in to single array (this seems to be what SK Hynix is suggesting) and then updated ind1 using the R value
// Fortran -> 1-indexed
// Ranged Read as I understand the SK Hynix implementation
// ind1_orig = [1, 2, 3, 6, 7, 8]
// buffer = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
// res_orig = [1.0, 2.0, 3.0, 6.0, 7.0, 8.0]

// R = 3
// ind1_transformed = [1, 6]
// buffer = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
// res_transformed = [1.0, 2.0, 3.0, 6.0, 7.0, 8.0]
// This requires merging/interleaving a few arrays in to a single buffer (buffer) before the call
// It also requires writing the loops to account for the merged/interleaved result (res)

// Example
// a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
// b = [11.0, 12.0, 13.0, 14.0, 15.0, 16.0]
// ind1_orig = [1, 2, 5, 6]
// res_a_packed = [1.0, 2.0, 5.0, 6.0]
// res_b_packed = [11.0, 12.0, 15.0, 16.0]

// Becomes
// a_b_merged = [1.0, 11.0, 2.0, 12.0, 3.0, 13.0, 4.0, 14.0, 5.0, 15.0, 6.0, 16.0]
// R = 2
// ind1 = [1, 3, 9, 11] -> ind1[i] = ind1_orig[i] * R - 1
// res_ab_merged_packed = [1.0, 11.0, 2.0, 12.0, 5.0, 15.0, 6.0, 16.0] -> length of the result is length(ind1) * 2 and can be calculated ahead of time
void scoria_read_ranged_1_c(double *res, const double *buffer, const size_t N, const size_t *ind1, const size_t N1, const size_t R) {
  size_t count = 0;
  for (size_t i = 0; i < N1; ++i) {
    for(size_t j = ind1[i] - 1; j < ind1[i] + R - 1; ++j, ++count) {
      res[count] = buffer[j];
    }
  }
}

void scoria_write_ranged_1_c(double *buffer, const double *input, const size_t N, const size_t *ind1, const size_t N1, const size_t R) {
  size_t count = 0;
  for (size_t i = 0; i < N1; ++i) {
    for (size_t j = ind1[i] - 1; j < ind1[i] + R - 1; ++j, ++count) {
      buffer[ind1[i] - 1] = input[count];
    }
  }
}