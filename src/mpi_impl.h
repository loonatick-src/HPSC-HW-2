#ifndef _MY_MPI_IMPL_H
#define _MY_MPI_IMPL_H
#include <stdio.h>

int
read_matrices(double *M_1, double *M_2,
        FILE *m1_file, FILE *m2_file, int width,
        int proc_rank, int num_procs);

int
matMulSquare_baseline_mpi(const double *M_1, const double *M_2,
        double *P, int width,
        int proc_rank, int num_procs);

int
matMulSquare_transpose_mpi(const double *M_1, const double *M_2,
        double *P, int width,
        int proc_rank, int num_procs);

int
matMulSquare_pretranspose_mpi(const double *M_1, const double *M_2,
        double *P, int width,
        int proc_rank, int num_procs);
#endif
