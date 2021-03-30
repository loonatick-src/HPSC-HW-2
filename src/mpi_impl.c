#include "dbg.h"
#include "mpi_impl.h"
#include <inttypes.h>
#include <stdlib.h>


int
matMulSquare_baseline(const double *M_1, const double *M_2, double *P, const uint32_t width)
{
    int mpi_init_flag;
    mpi_err = MPI_Init(NULL, NULL);
    check(mpi_err, "MPI environment failed to initialize");

    int proc_rank, num_procs;
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    check(mpi_err, "Call to `MPI_Comm_size` returned with error");
    check(num_procs > 2, "Mate, why even use MPI if your number of processes is %d?", num_procs);

    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    check(mpi_err, "Call to `MPI_Comm_rank` returned with error");

    const uint32_t matsize = width * width; 

     
    return EXIT_SUCCESS;

error:
    mpi_err = MPI_Initialized(&mpi_init_flag);

    if (mpi_err)
    {
        log_warn("Call to `MPI_Initialized` returned with error");
    }

    if (mpi_init_flag)
    {
        log_info("An MPI environment was initialized.\
                Attempting to finalize...");
        mpi_err = MPI_Finalize();
        if (mpi_err)
        {
            log_warn("MPI_Finalize() returned with error during cleanup");
        } else
        {
            log_info("MPI Finalized");
        }
    }

    return EXIT_FAILURE;
}


