#ifndef _MY_OMP_IMPL_H
#define _MY_OMP_IMPL_H
void
matMulSquare_baseline(const double *M_1,
                      const double *M_2,
                      double * P,
                      const uint32_t width);

int
matMulSquare_transpose(const double *M_1,
                       const double *M_2,
                       double * P,
                       const uint32_t width);
#endif
