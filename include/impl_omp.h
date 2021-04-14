#ifndef _MY_OMP_IMPL_H
#define _MY_OMP_IMPL_H

#include <stdint.h>

#define ARGUMENT_SIGNATURE_OMP const double *M_1, const double *M_2, double *P, uint32_t width


typedef int (*impl_omp_t)(ARGUMENT_SIGNATURE_OMP);


int
matMulSquare_baseline_omp(ARGUMENT_SIGNATURE_OMP);

int
matMulSquare_transpose_omp(ARGUMENT_SIGNATURE_OMP);

int
matMulSquare_pretranspose_omp(ARGUMENT_SIGNATURE_OMP);

int gaussian_elimination_naive_inplace_omp(double *M, uint32_t width);

#endif
