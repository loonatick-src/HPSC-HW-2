#include "dbg.h"
#include "impl_mpi.h"

#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 0.01;


int main(int argc, char *argv[])
{
    double *m1 = NULL, *m2 = NULL, *p = NULL;
    int width, proc_rank, num_procs;
    int mpi_err = MPI_Init(&argc, &argv);
    check(!mpi_err, "MPI failed to initialize.");

    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    check(!mpi_err, "MPI_Comm_size failed");
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    check(!mpi_err, "MPI_Comm_rank failed");

    if (proc_rank == 0)
    {
        int scan_count = scanf("%d", &width);
        check(scan_count != EOF, "Unexpected EOF");
        check(scan_count > 0, "Nothing scanned");
        int mat_size = width * width;

        m1 = (double *)malloc(mat_size * sizeof(double));
        m2 = (double *)malloc(mat_size * sizeof(double));
        p = (double *)malloc(mat_size * sizeof(double));
        for (int i = 0; i < mat_size; i++)
        {
            scan_count = scanf("%lf", m1 + i);
            check(scan_count != EOF, "Unexpected EOF");
            check(scan_count > 0, "Nothing scanned");
        }
        for (int i = 0; i < mat_size; i++)
        {
            scan_count = scanf("%lf", m2 + i);
            check(scan_count != EOF, "Unexpected EOF");
            check(scan_count > 0, "Nothing scanned");
        }
    }

    MPI_Bcast(&width, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    int mat_size = width * width;
    int my_err;

#ifdef BALANCED_BASE
    my_err = matMulSquare_balanced_mpi(m1, m2, p, width, proc_rank, num_procs);
#endif
#ifdef TRANSPOSE
    my_err = matMulSquare_transpose_mpi(m1, m2, p, width, proc_rank, num_procs);
#endif
#ifdef BASELINE
    my_err = matMulSquare_baseline_mpi(m1, m2, p, width, proc_rank, num_procs);
#endif
#ifdef PRETRANSPOSE
    my_err = matMulSquare_pretranspose_mpi(m1, m2, p, width, proc_rank, num_procs);
#endif
    check(!my_err, "Process %d: Something went wrong in matMul", proc_rank);
    if (proc_rank == 0) {
        for (int i = 0; i < mat_size; i++)
        {
            double elmt;
            scanf("%lf", &elmt);
            double residue = elmt - p[i];
            if (residue >= threshold)
            {
                log_warn("At index %d:\
                    expected %lf, eound %lf", i, elmt, p[i]);

            }
            else 
            {
                log_info("%d: Ye", i);
            }
            //(residue <= threshold, "At index %d:\
                    expected %lf, eound %lf", i, elmt, p[i]);
        }

        free(m1);
        free(m2);
        free(p);
    }
    debug_mpi(proc_rank, "Exit success");
    return EXIT_SUCCESS; 
error:
    if (m1) free(m1);
    if (m2) free(m2);
    if (p) free(p);
    return EXIT_FAILURE;
}

// mpi_bcast
