#include "dbg.h"
#include "impl_mpi.h"
#include "matrixio.h"

#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 0.01;

impl_mpi_t matmul_methods_mpi[] = {
                              matMulSquare_balanced_mpi,
                              matMulSquare_transpose_mpi,
                              matMulSquare_pretranspose_mpi};

const int num_methods_mpi = sizeof(matmul_methods_mpi)/sizeof(impl_mpi_t);

int main(int argc, char *argv[])
{
    double *m1 = NULL, *m2 = NULL, *p = NULL;
    int width, proc_rank, num_procs;
    int mpi_err, mpi_init_flag;
    int method_index = 0;  // defaults to baseline
    impl_mpi_t matmul;
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

    mpi_err = MPI_Bcast(&method_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    check(!mpi_err, "Call to MPI_Bcast returned with error");
    debug_mpi(proc_rank,"Method_index: %d", method_index);
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
        debug("matrix read complete");
    }

    mpi_err = MPI_Bcast(&width, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    check(!mpi_err, "`MPI_Bcast` returned with error.");
    int mat_size = width * width;

    debug_mpi(proc_rank,"Performing matmul");
    double start_time = MPI_Wtime();
    int my_err = matmul(m1, m2, p, width, proc_rank, num_procs);
    double end_time = MPI_Wtime();
    check(!my_err, "Something went wrong during matrix multiplication");
    debug_mpi(proc_rank, "Returned from matmul");

    double execution_time_matmul = end_time - start_time;

    if (proc_rank == 0) {
        for (int i = 0; i < mat_size; i++)
        {
            double elmt;
            int scan_rv = scanf("%lf", &elmt);
            check(scan_rv != EOF, "Unexpected EOF while reading matrix");
            check(scan_rv > 0, "Nothing was scanned");
            double residue = elmt - p[i];
            check(residue < threshold, "Bad numericals\
                    - matrix multiplication test case failed");
        }
        free(m1);
        free(m2);
    }
    debug_mpi(proc_rank, "m1 and m2 freed");

    start_time = MPI_Wtime();
    my_err = gaussian_elimination_naive_inplace_mpi(p, width, proc_rank, num_procs);
    end_time = MPI_Wtime();
    check(!my_err, "Error during gaussian elimination");
    debug("returned from gaussian elimination");
    double execution_time_elimination = end_time - start_time;

    //if (proc_rank == 0)
    //{
     //   for (int row = 0; row < width; row++)
      //  {
       //     for (int col = 0; col < width; col++)
        //    {
         //       printf("%lf ", p[row * width + col]);
          //  }
           // putchar('\n');
       // }
   //}

    // width num_procs matmul_time elim_time
    if (proc_rank == 0)
    {
        printf("%d %d %lf %lf\n", width, num_procs, 
                execution_time_matmul,
                execution_time_elimination); 
    }
    debug_mpi(proc_rank, "Exit success");
    if (p)
        free(p);
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
