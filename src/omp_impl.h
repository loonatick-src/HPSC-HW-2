#ifndef _MY_OMP_IMPL_H
#define _MY_OMP_IMPL_H
#include <stdint.h>
int
matMulSquare_baseline(const double *M_1,
                      const double *M_2,
                      double * P,
                      uint32_t width);

int
matMulSquare_transpose(const double *M_1,
                       const double *M_2,
                       double * P,
                       uint32_t width);

int
matMulSquare_pretranspose(const double *M_1,
                       const double *M_2,
                       double * P,
                       uint32_t width);

#endif
