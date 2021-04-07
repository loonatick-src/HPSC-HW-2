#include "dbg.h"
#include "impl_mpi.h"
#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NOT_IMPLEMENTED_ERROR 30


int
read_matrices(double *M_1, double *M_2,
        FILE *m1_file, FILE *m2_file, int width, 
        int proc_rank, int num_procs)
{
   /* reads matrices from the files into the malloc'd memory */
    const int matrix_size = width * width;
    int mpi_init_flag;
    int mpi_err = MPI_Initialized(&mpi_init_flag);
    check(!mpi_err, "Call to `MPI_Initialized` returned with error");
    debug("MPI_Initialized crossed");
    check(mpi_init_flag, "No MPI environment was found to be initialized");

    if (proc_rank == 0)
    {
    check_mem(M_1);
    check_mem(M_2);
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
    }

    return EXIT_SUCCESS;
error:
    // let caller handle error
    return EXIT_FAILURE;
}


int
matMulSquare_transpose_mpi(const double *M_1, double *M_2, double *P,
        int width, int proc_rank, int num_procs)
{
    const int num_rows_per_proc = width / num_procs;
    check(num_rows_per_proc > 0, 
            "Poorly balanced problem;\
            width of matrices: %d, number of processes: %d",
            width, num_procs);
    const int anomalous_num_rows = num_rows_per_proc + (width % num_procs);
    const int anomalous_proc = num_procs - 1;
    int num_elements_per_proc = num_rows_per_proc * width;
    const int anomalous_num_elements = anomalous_num_rows * width;
    if (proc_rank == anomalous_proc)
        num_elements_per_proc = anomalous_num_elements;
    const int mat_size = width * width;

    double *send_buffer = NULL;
    double *recv_buffer = NULL;
    double *M_2tr = (double *)malloc(mat_size * sizeof(double));
    check_mem(M_2tr);
    if (proc_rank == 0)
    {
        check_mem(M_1); check_mem(M_2); check_mem(P);
        for (int i = 0; i < mat_size; i++)
        {
            const int col = i%width;
            const int row = i/width;
            M_2tr[col*row + width] = M_2[row*col + width];
        }
    }

    // broadcasting the transposed matrix from process 0
    MPI_Bcast(M_2tr, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (proc_rank == 0)
    {
        int skip = 0;
        for (int proc = 1; proc < num_procs-1; proc++)
        {
            MPI_Send(M_1 + skip, num_elements_per_proc,
                    MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
            skip += num_elements_per_proc;
        }

        MPI_Send(M_1 + skip, anomalous_num_elements,
                MPI_DOUBLE, anomalous_proc, 0, MPI_COMM_WORLD);

        skip = 0;
        for (int proc = 1; proc < num_procs-1; proc++)
        {
            MPI_Recv(P + skip, num_elements_per_proc, MPI_DOUBLE,
                    proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            skip += num_elements_per_proc;
        }

        MPI_Recv(P + skip, anomalous_num_elements, MPI_DOUBLE,
                anomalous_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else 
    {
        recv_buffer = (double *)malloc(num_elements_per_proc);
        check_mem(recv_buffer);
        send_buffer = (double *)malloc(num_elements_per_proc);
        check_mem(send_buffer);
        MPI_Recv(recv_buffer, num_elements_per_proc, MPI_DOUBLE,
                0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < num_elements_per_proc; i++)
        {
            double sum = 0.0l;
            const int row = i/width;
            const int col = i%width;
            for (int j = 0; j < width; j++)
            {
                sum += recv_buffer[row*width + j] * M_2tr[col*width + j];
            }
            send_buffer[i] = sum;
        }
        MPI_Send(send_buffer, num_elements_per_proc, MPI_DOUBLE,
                0, 1, MPI_COMM_WORLD);
        free(recv_buffer);
        free(send_buffer);
    }

    free(M_2tr);
    return EXIT_SUCCESS;
error:
    if (M_2tr) free(M_2tr); 
    if (recv_buffer) free(recv_buffer);
    if (send_buffer) free(send_buffer);
    return EXIT_FAILURE;
}

int
matMulSquare_baseline_mpi(const double *M_1, double *M_2, double *P,
        int width, int proc_rank, int num_procs)
{
    const int num_rows_per_proc = width / num_procs;
    check(num_rows_per_proc > 0, 
            "Poorly balanced problem;\
            width of matrices: %d, number of processes: %d",
            width, num_procs);
    const int anomalous_num_rows = num_rows_per_proc + (width % num_procs);
    const int anomalous_proc = num_procs - 1;
    int num_elements_per_proc = num_rows_per_proc * width;
    const int anomalous_num_elements = anomalous_num_rows * width;
    if (proc_rank == anomalous_proc)
        num_elements_per_proc = anomalous_num_elements;
    const int mat_size = width * width;

    double *send_buffer = NULL;
    double *recv_buffer = NULL;
    if (proc_rank == 0)
    {
        check_mem(M_1); check_mem(M_2); check_mem(P);
        
    }

    // broadcasting the transposed matrix from process 0
    MPI_Bcast(M_2, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (proc_rank == 0)
    {
        int skip = 0;
        for (int proc = 1; proc < num_procs-1; proc++)
        {
            MPI_Send(M_1 + skip, num_elements_per_proc,
                    MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
            skip += num_elements_per_proc;
        }

        MPI_Send(M_1 + skip, anomalous_num_elements,
                MPI_DOUBLE, anomalous_proc, 0, MPI_COMM_WORLD);

        skip = 0;
        for (int proc = 1; proc < num_procs-1; proc++)
        {
            MPI_Recv(P + skip, num_elements_per_proc, MPI_DOUBLE,
                    proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            skip += num_elements_per_proc;
        }

        MPI_Recv(P + skip, anomalous_num_elements, MPI_DOUBLE,
                anomalous_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else 
    {
        recv_buffer = (double *)malloc(num_elements_per_proc);
        check_mem(recv_buffer);
        send_buffer = (double *)malloc(num_elements_per_proc);
        check_mem(send_buffer);
        MPI_Recv(recv_buffer, num_elements_per_proc, MPI_DOUBLE,
                0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < num_elements_per_proc; i++)
        {
            double sum = 0.0l;
            const int row = i/width;
            const int col = i%width;
            for (int j = 0; j < width; j++)
            {
                sum += recv_buffer[row*width + j] * M_2[j*width + col];
            }
            send_buffer[i] = sum;
        }
        MPI_Send(send_buffer, num_elements_per_proc, MPI_DOUBLE,
                0, 1, MPI_COMM_WORLD);
        free(recv_buffer);
        free(send_buffer);
    }

    return EXIT_SUCCESS;
error:
    if (recv_buffer) free(recv_buffer);
    if (send_buffer) free(send_buffer);
    return EXIT_FAILURE;
}


int 
matMulSquare_pretranspose_mpi(const double *M_1, double *M_2, double *P,
        int width, int proc_rank, int num_procs)
{
    const int num_rows_per_proc = width / num_procs;
    check(num_rows_per_proc > 0, 
            "Poorly balanced problem;\
            width of matrices: %d, number of processes: %d",
            width, num_procs);
    const int anomalous_num_rows = num_rows_per_proc + (width % num_procs);
    const int anomalous_proc = num_procs - 1;
    int num_elements_per_proc = num_rows_per_proc * width;
    const int anomalous_num_elements = anomalous_num_rows * width;
    if (proc_rank == anomalous_proc)
        num_elements_per_proc = anomalous_num_elements;
    const int mat_size = width * width;

    double *send_buffer = NULL;
    double *recv_buffer = NULL;
    if (proc_rank == 0)
    {
        check_mem(M_1); check_mem(M_2); check_mem(P);
        
    }

    // broadcasting the transposed matrix from process 0
    MPI_Bcast(M_2, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (proc_rank == 0)
    {
        int skip = 0;
        for (int proc = 1; proc < num_procs-1; proc++)
        {
            MPI_Send(M_1 + skip, num_elements_per_proc,
                    MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
            skip += num_elements_per_proc;
        }

        MPI_Send(M_1 + skip, anomalous_num_elements,
                MPI_DOUBLE, anomalous_proc, 0, MPI_COMM_WORLD);

        skip = 0;
        for (int proc = 1; proc < num_procs-1; proc++)
        {
            MPI_Recv(P + skip, num_elements_per_proc, MPI_DOUBLE,
                    proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            skip += num_elements_per_proc;
        }

        MPI_Recv(P + skip, anomalous_num_elements, MPI_DOUBLE,
                anomalous_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else 
    {
        recv_buffer = (double *)malloc(num_elements_per_proc);
        check_mem(recv_buffer);
        send_buffer = (double *)malloc(num_elements_per_proc);
        check_mem(send_buffer);
        MPI_Recv(recv_buffer, num_elements_per_proc, MPI_DOUBLE,
                0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0 ; i < num_elements_per_proc; i++)
        {
            double sum = 0.0l;
            const int row = i/width;
            const int col = i%width;
            for (int j = 0; j < width; j++)
            {
                sum += recv_buffer[row*width + j] * M_2[col*width + j];
            }
            send_buffer[i] = sum;
        }
        MPI_Send(send_buffer, num_elements_per_proc, MPI_DOUBLE,
                0, 1, MPI_COMM_WORLD);
        free(recv_buffer);
        free(send_buffer);
    }

    return EXIT_SUCCESS;
error:
    if (recv_buffer) free(recv_buffer);
    if (send_buffer) free(send_buffer);
    return EXIT_FAILURE;
}

