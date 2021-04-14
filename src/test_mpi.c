#include "dbg.h"
#include "impl_mpi.h"
#include "matrixio.h"

#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 0.01;

implementation_t matmul_methods_mpi[] = {matMulSquare_baseline_mpi,
                              matMulSquare_balanced_mpi,
                              matMulSquare_transpose_mpi,
                              matMulSquare_pretranspose_mpi};

const int num_methods_mpi = sizeof(matmul_methods_mpi)/sizeof(implementation_t);

int main(int argc, char *argv[])
{
    double *m1 = NULL, *m2 = NULL, *p = NULL;
    int width, proc_rank, num_procs;
    int mpi_err, mpi_init_flag;
    int method_index = 0;  // defaults to baseline
    implementation_t matmul;
    mpi_err = MPI_Init(&argc, &argv);
    check(!mpi_err, "MPI failed to initialize.");

    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    check(!mpi_err, "MPI_Comm_size failed");
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    check(!mpi_err, "MPI_Comm_rank failed");

    if (argc > 1)
    {
        method_index = strtol(argv[1], NULL, 10);
        if (method_index >= num_methods_mpi || method_index < 0)
        {
            log_warn("Invalid implementation number. Defaulting to 0");
            method_index = 0;
        }
    }
    matmul = matmul_methods_mpi[method_index];

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

    mpi_err = MPI_Bcast(&width, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    check(!mpi_err, "`MPI_Bcast` returned with error.");
    int mat_size = width * width;
    int my_err = matmul(m1, m2, p, width, proc_rank, num_procs);
    check(!my_err, "Something went wrong during matrix multiplication");

    if (proc_rank == 0) {
        for (int i = 0; i < mat_size; i++)
        {
            double elmt;
            scanf("%lf", &elmt);
            double residue = elmt - p[i];
            if (residue >= threshold)
            {
                log_warn("At index %d:\
                    expected %lf, found %lf", i, elmt, p[i]);

            }
            else 
            {
                log_info("%d: Ye", i);
            }
        }

        free(m1);
        free(m2);
    }

    my_err = gaussian_elimination_naive_inplace(p, width, proc_rank, num_procs);

    if (proc_rank == 0)
    {
        for (int row = 0; row < width; row++)
        {
            for (int col = 0; col < width; col++)
            {
                printf("%lf ", p[row * width + col]);
            }
            putchar('\n');
        }
    }
    

    debug_mpi(proc_rank, "Exit success");
    MPI_Finalize();
    return EXIT_SUCCESS; 
error:
    mpi_err = MPI_Initialized(&mpi_init_flag);
    if (mpi_err)
    {
        log_warn("Call to `MPI_Initialized` failed during cleanup");
    } else if (mpi_init_flag)
    {
        log_info("An MPI environment was found initialized. Finalizing");
        mpi_err = MPI_Finalize();
        if (mpi_err)
        {
            log_warn("Call to `MPI_Finalize` failed");
        }
    }
    if (m1) free(m1);
    if (m2) free(m2);
    if (p) free(p);
    return EXIT_FAILURE;
}

// mpi_bcast
