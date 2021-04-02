#include "dbg.h"
#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.b>
#include <string.h>


int
read_matrices(double *M_1, double *M_2, FILE *m1_file, FILE *m2_file, int width, int proc_rank, int num_procs);
{
    /* reads matrices from the files into the malloc'd memory */
    const int matrix_size = width * width;
    int mpi_init_flag;
    int mpi_err = MPI_Initialized(&mpi_init_flag);
    check(!mpi_err, "Call to `MPI_Initialized` returned with error");
    check(!mpi_init_flag, "No MPI environment was found to be initialized");

    check_mem(M_1);
    check_mem(M_2);
    if (proc_rank == 0)
    {
    check(m1_file, "Unopened file has been passed");
    check(m2_file, "Unopened file has been passed");

    for (int i = 0; i < matrix_size; i++)
    {
        int fscanf_ret = fscanf(m1_file, "%lf", &M_1[i]);
        check(fscanf_ret != EOF, "Unexpected EOF, left matrix. Scanned %d elements, promised %d elements", i, matrix_size);
        check(fscanf_ret != 0, "Failed to scan anything after scanning %d elements of the left matrix", i);
        fscanf_ret = fscanf(m2_file, "%lf", &M_2[i]);
        check(fscanf_ret != EOF, "Unexpected EOF, right matrix. Scanned %d elements, promised %d elements", i, matrix_size);
        check(fscanf_ret != 0, "Failed to scan anything after scanning %d elements of the right matrix", i);
    }

    return EXIT_SUCCESS;
error:
    mpi_err = MPI_Initialized(&mpi_init_flag);
    if (mpi_err)
    {
        log_warn("Call to `MPI_Finalized` returned with error");
    } else if (mpi_init_flag) 
    {
        log_info("An MPI environment was found to be initialized during cleanup. Attempting to finalize");
        mpi_err = MPI_Finalize();
        if (mpi_err)
            log_warn("Call to `MPI_Finalize` returned with error");
    }
    return EXIT_FAILURE;
}
