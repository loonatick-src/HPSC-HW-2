#include "dbg.h"
#include "impl_mpi.h"

#include <mpi.h>
#include <stdlib.h>

#define NOT_IMPLEMENTED 30

#ifdef DNDEBUG
#define debug_proc(R, M, ...) 
#else
#define debug_proc(R, M, ...) if (proc_rank == (R)) debug_mpi(R, M, ##__VA_ARGS__)
#endif

int
matMulSquare_baseline_mpi(const double *M_1, double *M_2,
        double *P, int width,
        int proc_rank, int num_procs)
{
    int num_rows_per_proc = width / num_procs;
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

    return EXIT_SUCCESS;
error:
    if (proc_rank != 0 && M_2)
        free(M_2);
    return EXIT_FAILURE;
}


int
matMulSquare_balanced_mpi(const double *M_1, double *M_2,
        double *P, int width,
        int proc_rank, int num_procs)
{
    const int mat_size = width * width;
    int num_rows_per_proc = width / num_procs;
    int unbalanced_num_rows = width - (num_rows_per_proc * num_procs);
    if (proc_rank == 0)
    {
        check_mem(M_1); check_mem(M_2); check_mem(P);
    }
    check(num_rows_per_proc > 0, "Poorly balanced problem: (%d rows, %d processes)", width, num_procs);

    int num_elements_per_proc = num_rows_per_proc * width;
    
    // all processes will need a copy of process M2
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

    // determine number of rows of M_1 to send to each process
    for (int i = 0; i < num_procs-1; i++)
    {
        send_counts[i] = num_elements_per_proc;
        if (unbalanced_num_rows != 0)
        {
            send_counts[i] += width;
            unbalanced_num_rows--;
        }
        displacements[i+1] = displacements[i] + send_counts[i];
    }
    send_counts[num_procs-1] = mat_size / num_procs;


    double *recv_buf = NULL, *send_buf = NULL;  
    recv_buf = (double *)malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(recv_buf);
    send_buf = (double *)malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(recv_buf);
     

    // scattering M_1
    MPI_Scatterv(M_1, send_counts, displacements, MPI_DOUBLE,
            recv_buf, send_counts[proc_rank], MPI_DOUBLE, 
            0, MPI_COMM_WORLD);


    num_rows_per_proc = send_counts[proc_rank]/width;
    for (int row = 0; row < num_rows_per_proc; row++)
    {
        for (int col = 0; col < width; col++)
        {
            int index = row*width + col;
            double sum = 0.0l;
            for (int i = 0; i < width; i++)
            {
                sum += recv_buf[row * width + i] * M_2[i * width + col];

            }
            send_buf[index] = sum;
        }
    }

    MPI_Gatherv(send_buf, send_counts[proc_rank], MPI_DOUBLE,
            P, send_counts, displacements, MPI_DOUBLE,
            0, MPI_COMM_WORLD);

    free(send_buf);
    free(recv_buf);

    
    if (proc_rank != 0)
    {
        free(M_2);
    }

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
    const int mat_size = width * width;
    int num_rows_per_proc = width / num_procs;
    int unbalanced_num_rows = width - (num_rows_per_proc * num_procs);

    if (proc_rank == 0)
    {
        check_mem(M_1);
        check_mem(M_2); 
        check_mem(P);
    }

    check(num_rows_per_proc > 0, "Poorly balanced problem: (%d rows, %d processes)", width, num_procs);

    int num_elements_per_proc = num_rows_per_proc * width;
    
    double *M_2tr = (double *)malloc(mat_size * sizeof(double)); 
    if (proc_rank == 0)
    {
        for (int i = 0; i < mat_size; i++)
        {
            const int row = i/width;
            const int col = i%width;
            M_2tr[col*width + row] = M_2[row*width + col];
        }
    }

    MPI_Bcast(M_2tr, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // send_counts for scatter_v
    int *send_counts = (int *)malloc(num_procs * sizeof(int));
    int *displacements = (int *)malloc(num_procs * sizeof(int));
    displacements[0] = 0;

    // determine number of rows of M_1 to send to each process
    for (int i = 0; i < num_procs-1; i++)
    {
        send_counts[i] = num_elements_per_proc;
        if (unbalanced_num_rows != 0)
        {
            send_counts[i] += width;
            unbalanced_num_rows--;
        }
        displacements[i+1] = displacements[i] + send_counts[i];
    }
    send_counts[num_procs-1] = mat_size / num_procs;


    double *recv_buf = NULL, *send_buf = NULL;  
    recv_buf = (double *)malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(recv_buf);
    send_buf = (double *)malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(recv_buf);
     

    // scattering M_1
    MPI_Scatterv(M_1, send_counts, displacements, MPI_DOUBLE,
            recv_buf, send_counts[proc_rank], MPI_DOUBLE, 
            0, MPI_COMM_WORLD);


    num_rows_per_proc = send_counts[proc_rank]/width;
    for (int row = 0; row < num_rows_per_proc; row++)
    {
        for (int col = 0; col < width; col++)
        {
            int index = row*width + col;
            double sum = 0.0l;
            for (int i = 0; i < width; i++)
            {
                sum += recv_buf[row * width + i] * M_2tr[col * width + i];

            }
            send_buf[index] = sum;
        }
    }

    MPI_Gatherv(send_buf, send_counts[proc_rank], MPI_DOUBLE,
            P, send_counts, displacements, MPI_DOUBLE,
            0, MPI_COMM_WORLD);

    free(M_2tr);
    free(send_buf);
    free(recv_buf);

    return EXIT_SUCCESS;
error:
    if (proc_rank != 0 && M_2)
        free(M_2);
    return EXIT_FAILURE;
}


int
matMulSquare_pretranspose_mpi(const double *M_1, double *M_2,
        double *P, int width,
        int proc_rank, int num_procs)
{
    const int mat_size = width * width;
    int num_rows_per_proc = width / num_procs;
    int unbalanced_num_rows = width - (num_rows_per_proc * num_procs);

    if (proc_rank == 0)
    {
        check_mem(M_1);
        check_mem(M_2); 
        check_mem(P);
    }

    check(num_rows_per_proc > 0, "Poorly balanced problem: (%d rows, %d processes)", width, num_procs);

    int num_elements_per_proc = num_rows_per_proc * width;
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

    // determine number of rows of M_1 to send to each process
    for (int i = 0; i < num_procs-1; i++)
    {
        send_counts[i] = num_elements_per_proc;
        if (unbalanced_num_rows != 0)
        {
            send_counts[i] += width;
            unbalanced_num_rows--;
        }
        displacements[i+1] = displacements[i] + send_counts[i];
    }
    send_counts[num_procs-1] = mat_size / num_procs;


    double *recv_buf = NULL, *send_buf = NULL;  
    recv_buf = (double *)malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(recv_buf);
    send_buf = (double *)malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(recv_buf);
     

    // scattering M_1
    MPI_Scatterv(M_1, send_counts, displacements, MPI_DOUBLE,
            recv_buf, send_counts[proc_rank], MPI_DOUBLE, 
            0, MPI_COMM_WORLD);


    num_rows_per_proc = send_counts[proc_rank]/width;
    for (int row = 0; row < num_rows_per_proc; row++)
    {
        for (int col = 0; col < width; col++)
        {
            int index = row*width + col;
            double sum = 0.0l;
            for (int i = 0; i < width; i++)
            {
                sum += recv_buf[row * width + i] * M_2[col * width + i];

            }
            send_buf[index] = sum;
        }
    }

    MPI_Gatherv(send_buf, send_counts[proc_rank], MPI_DOUBLE,
            P, send_counts, displacements, MPI_DOUBLE,
            0, MPI_COMM_WORLD);

    if (proc_rank != 0)
        free(M_2);
    free(send_buf);
    free(recv_buf);

    return EXIT_SUCCESS;
error:
    if (send_buf)
        free(send_buf);
    if (recv_buf)
        free(recv_buf);
    if (proc_rank != 0 && M_2)
        free(M_2);
    return EXIT_FAILURE;
}


int
gaussian_elimination_naive_inplace(double *M, /*double *P,*/ int width,
        int proc_rank, int num_procs)
{
    if (proc_rank == 0)
        check_mem(M);
    double *pivot_buf = NULL, *proc_buf = NULL;
    int *send_counts = NULL, *displacements = NULL;

    send_counts = (int *) malloc(sizeof(int) * num_procs);
    check_mem(send_counts);
    displacements = (int *) malloc(sizeof(int) * num_procs);
    check_mem(send_counts);

    int rows_per_proc = width / num_procs;
    int unbalanced_num_rows = width - num_procs * (width / num_procs);
    displacements[0] = 0;

    for (int i = 0; i < num_procs-1; i++)
    {
        send_counts[i] = rows_per_proc * width;
        if (unbalanced_num_rows != 0)
        {
            send_counts[i] += width;
            unbalanced_num_rows--;
        }
        displacements[i+1] = displacements[i] + send_counts[i];
    }
    send_counts[num_procs-1] = rows_per_proc * width;

    pivot_buf = (double *) malloc(width * sizeof(double));
    check_mem(pivot_buf);
    proc_buf = (double *) malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(proc_buf);

    // scattering the matrix to all processes
    int mpi_err = MPI_Scatterv(M, send_counts, displacements,
            MPI_DOUBLE, proc_buf, send_counts[proc_rank],
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    check(!mpi_err, "Call to `MPI_Scatterv` returned with error");

    // iterating over rows of the matrix
    // each row except the last becomes the pivot row
    int pivot_row = 0, pivot_proc = 0;
    double *temp = NULL;
    while (pivot_row < width - 1)
    {
        if (proc_rank == pivot_proc)
        {
            temp = pivot_buf;
            pivot_buf = proc_buf + pivot_row * width - displacements[proc_rank];
            check_mem(pivot_buf);
            check_mem(temp);
        }

        MPI_Bcast(pivot_buf, width, MPI_DOUBLE,
                pivot_proc, MPI_COMM_WORLD);

        int num_rows = send_counts[proc_rank]/width;
        for (int row = 0; row < num_rows; row++)
        {
            int global_row = displacements[proc_rank]/width + row;
            if (global_row > pivot_row)
            {
                check_mem(pivot_buf);
                double pivot = pivot_buf[pivot_row];
                check(pivot != 0, "Singular pivot");
                double leverage = proc_buf[width*row + pivot_row];
                proc_buf[width * row + pivot_row] = 0.0l;
                double factor = leverage / pivot;
                for (int col = pivot_row + 1; col < width; col++)
                {
                    int index = row * width + col;
                    proc_buf[index] -=  pivot_buf[index] * factor;
                }
            }
        }
        pivot_row++; 
        if (proc_rank == pivot_proc)
        {
            pivot_buf = temp;
            check_mem(pivot_buf);
            temp = NULL;
        }

        if (pivot_proc < num_procs-1) 
        {
            if ((pivot_row * width)%(displacements[pivot_proc+1]) ==  0)
            {
                pivot_proc++;
            }
        }
    }

    mpi_err = MPI_Gatherv(proc_buf, send_counts[proc_rank],
            MPI_DOUBLE, M, send_counts,
            displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    check(!mpi_err, "Gathering to root failed");

    free(pivot_buf);
    free(proc_buf);
    free(displacements);
    free(send_counts);
    return 0;
error:
    if (pivot_buf)
        free(pivot_buf);
    if (proc_buf)
        free(proc_buf);
    if (send_counts)
        free(send_counts);
    if (displacements)
        free(displacements);
    return -1;
}
