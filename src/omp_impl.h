#ifndef _MY_OMP_IMPL_H
#define _MY_OMP_IMPL_H
void
matMulSquare_baseline(const double * const M_1,
                      const double * const M_2,
                      double * P,
                      const uint32_t width);

int
matMulSquare_transpose(const double * const M_1,
                       const double * const M_2,
                       double * P,
                       const uint32_t width);
#endif
