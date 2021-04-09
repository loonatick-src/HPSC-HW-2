#include "dbg.h"
#include "impl_mpi.h"

#include <mpi.h>
#include <stdlib.h>

#define NOT_IMPLEMENTED 30

int
read_matrices(double *M_1, double *M_2,
        FILE *m1_file, FILE *m2_file, int width,
        int proc_rank, int num_procs)
{
    return NOT_IMPLEMENTED;
}

int
matMulSquare_baseline_mpi(const double *M_1, double *M_2,
        double *P, int width,
        int proc_rank, int num_procs)
{
    int num_rows_per_proc = width / num_procs;
    const int unbalanced_proc = num_procs - 1;
    if (proc_rank == 0)
    {
        check_mem(M_1); check_mem(M_2); check_mem(P);
    }
    check(num_rows_per_proc > 0, "Poorly balanced problem: (%d rows, %d processes)", width, num_procs);
    int num_elements_per_proc = num_rows_per_proc * width;
    const int mat_size = width * width;
    // not worrying about load balancing for the time being
    const int unbalanced_num_elements = mat_size - ((num_elements_per_proc) * (num_procs-1));
    const int unbalanced_num_rows = unbalanced_num_elements/width;
    

    
    // all processes will need a copy of process 2
    if (proc_rank != 0)
    {
        M_2 = (double *)malloc(mat_size * sizeof(double));
        check_mem(M_2);
    }

    MPI_Bcast(M_2, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // send_counts for scatter_v
    int *send_counts = (int *)malloc(num_procs * sizeof(int));
    int *displacements = (int *)malloc(num_procs * sizeof(int));
    displacements[0] = 0;

    for (int i = 0; i < num_procs-1; i++)
    {
        send_counts[i] = num_elements_per_proc;
        displacements[i+1] = displacements[i] + num_elements_per_proc;
    }
    // last process gets more elements than the others
    send_counts[num_procs-1] = unbalanced_num_elements;

    if (proc_rank == num_procs - 1) 
    {
        num_elements_per_proc = unbalanced_num_elements;
        num_rows_per_proc = unbalanced_num_rows;
    }

    double *recv_buf = NULL, *send_buf = NULL;  
    recv_buf = (double *)malloc(num_elements_per_proc * sizeof(double));
    check_mem(recv_buf);
    send_buf = (double *)malloc(num_elements_per_proc * sizeof(double));
    check_mem(recv_buf);
     

    // scattering M_1
    MPI_Scatterv(M_1, send_counts, displacements, MPI_DOUBLE,
            recv_buf, num_elements_per_proc, MPI_DOUBLE, 
            0, MPI_COMM_WORLD);


    for (int row = 0; row < num_rows_per_proc; row++)
    {
        for (int col = 0; col < width; col++)
        {
            int index = row*width + col;
            double sum = 0.0l;
            for (int i = 0; i < width; i++)
            {
                sum += recv_buf[row * width + i] * M_2[i * width + col];
                if (proc_rank == unbalanced_proc)
                {
                    debug("%lf", sum);
                }

            }
            send_buf[index] = sum;
        }
    }

    MPI_Gatherv(send_buf, num_elements_per_proc, MPI_DOUBLE,
            P, send_counts, displacements, MPI_DOUBLE,
            0, MPI_COMM_WORLD);

    free(send_buf);
    free(recv_buf);

    
    if (proc_rank != 0)
    {
        free(M_2);
    }
    MPI_Finalize();

    return EXIT_SUCCESS;
error:
    if (proc_rank != 0 && M_2)
        free(M_2);
    return EXIT_FAILURE;
}

int
matMulSquare_transpose_mpi(const double *M_1, double *M_2,
        double *P, int width,
        int proc_rank, int num_procs)
{
    return NOT_IMPLEMENTED;
}

int
matMulSquare_pretranspose_mpi(const double *M_1, double *M_2,
        double *P, int width,
        int proc_rank, int num_procs)
{
    return NOT_IMPLEMENTED;
}

