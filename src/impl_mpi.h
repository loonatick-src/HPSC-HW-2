#ifndef _MY_MPI_IMPL_H
#define _MY_MPI_IMPL_H
#include <stdio.h>

#define ARGUMENT_SIGNATURE const double *M_1, double *M_2,\
    double *P, int width,\
    int proc_rank, int num_procs

int
matMulSquare_baseline_mpi(ARGUMENT_SIGNATURE);

int
matMulSquare_transpose_mpi(ARGUMENT_SIGNATURE);

int
matMulSquare_pretranspose_mpi(ARGUMENT_SIGNATURE);

int
matMulSquare_balanced_mpi(ARGUMENT_SIGNATURE);
#endif
