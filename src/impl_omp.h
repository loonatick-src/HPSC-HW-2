#ifndef _MY_OMP_IMPL_H
#define _MY_OMP_IMPL_H
#include <stdint.h>
int
matMulSquare_baseline_omp(const double *M_1,
                      const double *M_2,
                      double * P,
                      uint32_t width);

int
matMulSquare_transpose_omp(const double *M_1,
                       const double *M_2,
                       double * P,
                       uint32_t width);

int
matMulSquare_pretranspose_omp(const double *M_1,
                       const double *M_2,
                       double * P,
                       uint32_t width);

#endif
