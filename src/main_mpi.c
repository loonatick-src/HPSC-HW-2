#include "dbg.h"
#include "mpi_impl.h"

#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define UINT64_TO_INT32_MASK 0x00000000ffffffff
#define INT32_SIGN_MASK 0x80000000

/* TODO
 * use getopt to get provide user with the option to input
   transposed right matrix for fast multiplication
 * complete the rest of the code
 */


int main(int argc, char *argv[])
{
    double *M_1, *M_2, *P;
    FILE *m1_file, *m2_file;
    int mpi_init_flag, mpi_err = MPI_Init(&argc, &argv);
    check(mpi_err, "Failed to initialize MPI");

    int proc_rank, num_procs;
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    check(mpi_err, "`Call to MPI_Comm_size` returned with error");
    check(num_procs != 1, "No, no base case for you.");
    check(num_procs > 2, "Mate, why even use MPI if you are not going to spawn more than 2 processes?");

    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    check(mpi_err, "Call to `MPI_Comm_rank` returned with error");

    if (proc_rank == 0) 
    {
        check(argc >= 4, "Usage: ./%s <n> <m1_path> <m2_path>,\
                where\nn: width of matrices (base 10)\
                \nm1_path: path of matrix 1 file (row-major)\
                \nm2_path: path of matrix 2 file (row-major)", argv[0]);

        const uint64_t widthl = strtoul(argv[1], NULL, 10);
        const uint64_t width_proto = UINT64_TO_INT32_MASK & widthl;
        check(width_proto == widthl, "Integer overflow");
        check(width_proto & INT32_SIGN_MASK > 0, "Integer overflow");

        const int width = (int) widthl;
        const int matsize = widthl * widthl;
        check(matsize >= widthl, "Integer overflow (number of matrix elements) - proposed matrix sizes are too largs");

        const int mem_size = matsize * sizeof(double);
        check(mem_size > matsize, "Integer overflow (number of bytes) - proposed matrix sizes are too large");

        const char *m1_path = argv[2];
        const char *m2_path = argv[3];

        m1_file = fopen(m1_path, "r");
        check(m1_file, "Failed to open %s", m1_path);

        m2_file = fopen(m2_path, "r");
        check(m2_file, "Failed to open %s", m2_path);

        M_1 = (double *)malloc(mem_size);
        check_mem(M_1);

        M_2 = (double *)malloc(mem_size);
        check_mem(M_2);

        P = (double *)malloc(mem_size);
        check_mem(P);

        int in_count = 0;
        double melmnt;
        while (in_count < matsize)
        {
            int fscanf_err = fscanf(m1_file, "%lf", &melmnt);
            check(scanf_err != EOF, "Unexpected EOF: scanned %lu elements, promised %lu elements", in_count, matsize);
            check(scanf_err, "`scanf` returned with error after scanning %lu elements", in_count);

            M_1[in_count++] = melmnt;
        }
        fclose(m1_file);

        in_count = 0;
        while (in_count < matsize)
        {
            int fscanf_err = fscanf(m2_file, "%lf", &melmnt);
            check(scanf_err != EOF, "Unexpected EOF: scanned %lu elements, promised %lu elements", in_count, matsize);
            check(scanf_err, "`scanf` returned with error after scanning %lu elements", in_count);

            M_2[in_count++] = melmnt;
        }
        fclose(m2_file);

        const int num_rows_per_proc = width / num_procs;
        check(num_rows_per_proc > 0, "Poorly balanced inputs.\
                \nWidth: %d, Number of MPI Processes: %d", width, num_procs);
        const int anomalous_proc = num_procs - 1;
        const int anomalous_num_rows = num_rows_per_proc + (num_rows % num_rows_per_proc);

        if ((num_procs-1) * num_rows_per_proc + anomalous_num_rows != width)
        {
            sentinel("There is a mistake in the logic\
                    for calculating number of rows to be sent to each process");
        }
        const int anomalous_num_elmts = anomalous_num_rows * width;


        const int num_elmts = width * num_rows_per_proc;
        for (int proc = 1; proc < num_procs-1; proc++)
        {
            MPI_Send(M_1 + (proc-1)*num_elmts,
                    num_elmts,
                    MPI_DOUBLE,
                    proc, 0,
                    MPI_COMM_WORLD);
        }
        MPI_Send(M_1 + anomalous_proc * num_elmts,
                anomalous_num_elmts, 
                MPI_DOUBLE,
                anomalous_proc, 0, 
                MPI_COMM_WORLD);
    } else
    {
        
    }


    free(M_1);
    free(M_2);
    free(P);

    mpi_err = MPI_Finlize();
    check(mpi_err, "MPI_Failed to finalize after everything else is done and dusted");
    return EXIT_SUCCESS;
error:
    

    if (M_1)
    {
        log_info("Freeing memory allocated for left input matrix\
                in process %d", proc_rank);
        free(M_1);
    }

    if (M_2)
    {
        log_info("Freeing memory allocated for right input matrix\
                in process %d", proc_rank);
    }

    if (P)
    {
        log_info("Freeing memory allocated for output matrix\
                in process %d", proc_rank);
    }

    if (m1_file)
    {
        fclose(m1_file);
    }

    if (m2_file)
    {
        fclose(m2_file);
    }

    mpi_err = MPI_Initialized(&mpi_init_flag);
    if (mpi_err)
    {
        log_warn("Call to `MPI_Initialized` returned error\
                during cleanup");
    } else if (mpi_init_flag)
    {
        log_info("MPI environment found running during cleanup.\
                attempting to finalize...");
        mpi_err = MPI_Finalize();
        if (mpi_err)
        {
            log_warn("MPI Environment failed to finalize");
        } else
        {
            log_info("MPI Environment finalized");
        }
    }    
    return EXIT_FAILURE;
}
