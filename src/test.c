#include "dbg.h"

#ifdef _MPI_TEST
#include <mpi.h>
#include "mpi_impl.h"
#endif
#ifdef _OMP_TEST
#include <omp.h>
#include "omp_impl.h"
#endif
#include <math.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE -1
#endif

const uint32_t width = 100;
const uint32_t matsize = width * width;
const uint32_t memsize = matsize * sizeof(double);
const double threshold_percent_error = 1.e-2;
const char *m1_path = "m_1_100.dat";
const char *m2_path = "m_2_100.dat";
const char *mmulpath = "matmul_100.dat";
const char *mmulpatht = "matmul_100t.dat";


int
read_matrix(double *m, FILE *filehandle, uint32_t matsize)
{
    check(matsize > 0, "Matrix size must be a positive integer, duh");
    for (int i = 0; i < matsize; i++)
    {
        int fscanf_out = fscanf(filehandle, "%lf", &(m[i]));
        check(fscanf_out != EOF, "There is nothing more to read from stdin");
        check(fscanf_out, "fscanf couldn't read anything");
    }
    return EXIT_SUCCESS;
error:
    return EXIT_FAILURE;
    
}


double
percent_error(double of, double against)
{
    debug("%lf %lf", of, against);
    return (fabs(of - against)/against) * 100.0l;
}


int
validate_matrix(double *v, double *a, uint32_t matsize)
{
    for (uint32_t index = 0; index < matsize; index++)
    {
        double err = percent_error(v[index], a[index]);
        if (err > threshold_percent_error)
        {
            log_err("Erroneous value at index %d, expected close to %lf, found %lf", index, a[index], v[index]);
            return -1;
        }
        index++;
    }
    return 0;
}


#ifdef _OMP_TEST
void
test_matmul_omp()
{
#ifdef _OPENMP
    omp_set_num_threads(6);
#endif
    log_info("Testing OpenMP matrix multiplication implementations");
    log_info("Allocating space for matrices");
    double *m1 = (double *)malloc(memsize);
    double *m2 = (double *)malloc(memsize);
    double *matmul = (double *)malloc(memsize);
    double *p = (double *)malloc(memsize);

    FILE *m1_file = NULL;
    FILE *m2_file = NULL;
    FILE *matmul_file = NULL;

    log_info("Opening files for reading matrices");
    m1_file = fopen(m1_path, "r");
    m2_file = fopen(m2_path, "r");
    matmul_file = fopen(mmulpath, "r");

    log_info("Reading matrices into memory");
    read_matrix(m1, m1_file, matsize);
    fclose(m1_file);

    read_matrix(m2, m2_file, matsize);
    fclose(m2_file);

    read_matrix(matmul, matmul_file, matsize);
    fclose(matmul_file);

    // First implementation
    log_info("Testing first implementation (direct multiplication)");
    int my_err = matMulSquare_baseline(m1, m2, p, width);
    check(!my_err, "Something went wrong in the first implementation");
    my_err = validate_matrix(p, matmul, matsize);
    check_debug(!my_err, "Output of first implementation seems erroneous");
    log_info("Testing of first implementation complete");


    // Second implementation
    log_info("Testing second implementation (transposed multiplication)");
    my_err = matMulSquare_transpose(m1, m2, p, width);
    check(!my_err, "Something went wrong in the second implementation");
    my_err = validate_matrix(p, matmul, matsize);
    check(!my_err, "Output of second implementation seems erroneous");
    log_info("Testing of second implementation complete");

    // Third implementation
    log_info("Reading transposed matrices into memory");
    matmul_file = fopen(mmulpatht, "r");
    check(matmul_file, "Error while opening %s: file may not exist", mmulpatht);

    my_err = read_matrix(matmul, matmul_file, matsize);
    check(!my_err, "Something went wrong while reading transposed product matrix");
    fclose(matmul_file);

    log_info("Testing third implementation (pretransposed multiplication");
    my_err = matMulSquare_pretranspose(m1, m2, p, width);
    check(!my_err, "Something went wrong with the third implementation");
    my_err = validate_matrix(p, matmul, matsize);
    check(!my_err, "Third implementation seems erroneous");
    log_info("Third implementation tested successfully");

#ifdef _OPENMP
    log_info("Testing of OpenMP implementations completed successfully");
#endif

    free(m1);
    free(m2);
    free(p);
    free(matmul);
    return ;
error:
    if (m1) free(m1);
    if (m2) free(m2);
    if (p) free(p);
    if (m1_file) fclose(m1_file);
    if (m2_file) fclose(m2_file);
}
#endif


#ifdef _MPI_TEST
void
test_matmul_mpi()
{
    // TODO: change mpi_err to mpi_code
    double *M_1 = NULL, *M_2 = NULL, *P = NULL;
    FILE *m1_file = NULL, *m2_file = NULL;

    int mpi_init_flag, mpi_err = MPI_Init(NULL, NULL);
    check(!mpi_err, "MPI Failed to initialize");

    int num_procs = 0, proc_rank = 0;
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    check(num_procs > 0, "Call to MPI_Comm_size failed");
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    check(num_procs > 0, "Call to MPI_Comm_rank failed");
    debug("%d: num_procs = %d", proc_rank, num_procs); 
    const int n = 100;
    const int nn = n * n;
    if (proc_rank == 0) 
    {
    M_1 = (double *) malloc(nn * sizeof(double));
    check_mem(M_1);
    M_2 = (double *) malloc(nn * sizeof(double));
    check_mem(M_2);
    P = (double *) malloc(nn * sizeof(double));
    check_mem(P);
    m1_file = fopen(m1_path, "r");
    m2_file = fopen(m2_path, "r");
    }

    log_info("Process %d: Testing `read_matrices`", proc_rank);
    int my_err = read_matrices(M_1, M_2, m1_file, m2_file, n,
            proc_rank, num_procs);
    check(!my_err, "Call to `read_matrices` returned with error");
    

    debug("%d: Testing baseline implementation", proc_rank);
    debug("%d: Freeing resources", proc_rank);
    if (M_1) free(M_1);
    if (M_2) free(M_2);
    if (P) free(P);
    if (m1_file) fclose(m1_file);
    if (m2_file) fclose(m2_file);
    mpi_err = MPI_Finalize();
    return;
error:
    if (m1_file)
        fclose(m1_file);
    if (m2_file)
        fclose(m2_file);
    if (M_1)
        free(M_1);
    if(M_2)
        free(M_2);
    if (P)
        free(P);

    mpi_err = MPI_Initialized(&mpi_init_flag);
    if (mpi_err)
    {
        log_warn("Call to `MPI_Initialized` failed");
    } else if (mpi_init_flag)
    {
        log_info("An MPI environment was found initialized during cleanup. Attempting to finalize...");
        mpi_err = MPI_Finalize();
        if (mpi_err)
        {
            log_warn("Call to `MPI_Finalize` failed");
        }
    }
}
#endif


int main(int argc, char *argv[])
{
#ifdef _OMP_TEST
    test_matmul_omp();
#endif

#ifdef _MPI_TEST
    test_matmul_mpi();
#endif
    return 0;
}
