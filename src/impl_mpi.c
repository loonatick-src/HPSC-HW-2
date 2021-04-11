#include "dbg.h"
#include "impl_mpi.h"

#include <mpi.h>
#include <stdlib.h>

#define NOT_IMPLEMENTED 30


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
        if (proc_rank == 0)
        {
            debug("send_counts[%d] = %d, displacements[%d] = %d", i, send_counts[i], i, displacements[i]);
        }

    }
    send_counts[num_procs-1] = mat_size / num_procs;
    if (proc_rank == 0)
    {
        debug("send_counts[%d] = %d", num_procs - 1, send_counts[proc_rank]);
    }


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
gaussian_elimination_naive(double *M, /*double *P,*/ int width,
        int proc_rank, int num_procs)
{
    int mpi_init_flag;
    int *send_counts = NULL, *displacements = NULL;
    double *send_buf = NULL, *recv_buf = NULL;
    int mpi_err = MPI_Initialized(&mpi_init_flag);
    check(!mpi_err, "`call to `MPI_Initialized` returned with error. Aborting.");
    check(mpi_init_flag, "MPI not initialized. Aborting.");
    check(width > 0, "Non-positive width.");

    if (proc_rank == 0)
    {
        check_mem(M);
        // check_mem(P);
    }

    send_counts = (int *)malloc(sizeof(int) * num_procs);
    displacements = (int *)malloc(sizeof(int) * num_procs);
    displacements[0] = 0;

    int base_num_rows = width / num_procs;
    int unbalanced_num_rows = width % num_procs;

    for (int proc = 0; proc < num_procs-1; proc++)
    {
        send_counts[proc] = base_num_rows * width;
        if (unbalanced_num_rows != 0)
        {
            send_counts[proc] += width;
            unbalanced_num_rows--;
        }
        displacements[proc+1] = displacements[proc] + send_counts[proc];
    }
    send_counts[num_procs-1] = base_num_rows * width;


    recv_buf = (double *)malloc(send_counts[proc_rank] * sizeof(double));
    check_mem(recv_buf);
        

    mpi_err = MPI_Scatterv(M, send_counts, displacements,
            MPI_DOUBLE, recv_buf, send_counts[proc_rank],
            MPI_DOUBLE, 0, MPI_COMM_WORLD);  

    //for (int i = 0; i < send_counts[proc_rank]; i++)
    //{
    //   debug_mpi(proc_rank, "recv_buf[%d] = %lf", i, recv_buf[i]);
    //}

    
    send_buf = (double *)malloc(width * sizeof(double));
    int iter = 0;
    int proc = 0;
    while (proc < num_procs);
    {
        for (int row = 0; row < send_counts[proc]/width; row++)
        {
            void *temp = NULL; 
            if (proc_rank == proc)
            {
                temp = send_buf;
                send_buf = recv_buf + width*row;
            }

            // broadcasting the pivot row
            MPI_Bcast(send_buf, width, MPI_DOUBLE,
                    proc, MPI_COMM_WORLD);

            double pivot = send_buf[iter];
            // TODO: do threshold check instead
            check(pivot != 0.0l, "Zero pivot found. Aborting (use partial pivoting algorithm).")
            if (proc_rank > proc)
            {
                for (int i = 0; i < send_counts[proc_rank]/width; i++)
                {
                    double factor = recv_buf[i * width + iter]/pivot;
                    recv_buf[i*width + iter] = 0.0l;
                    for (int j = iter + 1; j < width; j++)
                    {
                        recv_buf[i * width + j] -= send_buf[j]*factor;
                    }
                }
            }
            if (proc_rank == proc)
            {
                send_buf = temp;
            }
            iter++;
            if (iter == width-1)
                break;
        }
        proc++;
    }

    free(send_counts);
    free(displacements);
    free(recv_buf);
    return 0;
error:
    if (send_counts)
        free(send_counts);
    if (displacements)
        free(displacements);
    if (send_buf)
        free(send_buf);
    if (recv_buf)
        free(recv_buf);
    return -1;
}
