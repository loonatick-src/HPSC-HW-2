#include "dbg.h"
#include "impl_mpi.h"
#include "impl_omp.h"

#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

const double THRESHOLD = 0.01l;

static
inline
double percent_error(double of, double against)
{
    if (against != 0)
        return fabs(against - of)/fabs(against);
    else if (against == of)
        return 0.0l;
    else
        return -1.0l;
}


impl_mpi_t mpi_methods[] = 
{ 
    matMulSquare_balanced_mpi,
    matMulSquare_transpose_mpi,
    matMulSquare_pretranspose_mpi
};

impl_omp_t omp_methods[] = 
{
    matMulSquare_baseline_omp,
    matMulSquare_transpose_omp,
    matMulSquare_pretranspose_omp
};

const int num_methods = sizeof(omp_methods)/sizeof(impl_omp_t);

int main(int argc, char *argv[])
{
    double *m1 = NULL, *m2 = NULL;
    double *p_omp = NULL, *p_mpi = NULL;
    int mpi_err, scan_rv, my_err, mpi_init_flag;
    uint32_t width_omp; int width_mpi;

    mpi_err = MPI_Init(&argc, &argv);
    check(!mpi_err, "MPI initialilzation failed");

    int method_index = 0;  // defaults to baseline
    int proc_rank, num_procs;
    //int num_threads = omp_get_max_threads();

    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    check(!mpi_err, "MPI_Comm_size returned with error");
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    check(!mpi_err, "MPI_Comm_rank returned with error");
    
    if (proc_rank == 0)
    {
        if (argc > 1)
        {
            method_index = (int) strtol(argv[1], NULL, 10);
            if (method_index < 0 || method_index >= num_methods)
            {
                log_warn("Invalid method index %d, defaulting to 0.", method_index);
                method_index = 0;
            }
        }

        scan_rv = scanf("%d", &width_mpi);
        check(scan_rv != EOF, "Unexpected EOF while scanning width");
        check(scan_rv > 0, "Couldn't read a valid unsigned int value");
        check(width_mpi > 0, "Negative width/overflow");

        const int mat_size = width_mpi * width_mpi;
        m1 = (double *) malloc(mat_size * sizeof(double));
        check_mem(m1);
        m2 = (double *) malloc(mat_size * sizeof(double));
        check_mem(m2);
        p_mpi = (double *) malloc(mat_size * sizeof(double));
        check_mem(p_mpi);
        p_omp = (double *) malloc(mat_size * sizeof(double));
        check_mem(p_omp);

        for (int i = 0; i < mat_size; i++)
        {
            scan_rv = scanf("%lf", m1 + i);
            check(scan_rv != EOF, "Unexpected EOF");
            check(scan_rv > 0, "Could not scan a valid double");
        }
        for (int i = 0; i < mat_size; i++)
        {
            scan_rv = scanf("%lf", m2 + i);
            check(scan_rv != EOF, "Unexpected EOF");
            check(scan_rv > 0, "Could not scan a valid double");
        }
    }

    mpi_err = MPI_Bcast(&width_mpi, 1, MPI_INTEGER,
            0, MPI_COMM_WORLD);
    check(!mpi_err, "Broadcasting width failed");
    mpi_err = MPI_Bcast(&method_index, 1, MPI_INTEGER,
            0, MPI_COMM_WORLD);
    check(!mpi_err, "Broadcasting method_index failed");

    width_omp = (uint32_t) width_mpi;
    size_t width = (size_t) width_omp;
    impl_omp_t matmul_omp = omp_methods[method_index];
    impl_mpi_t matmul_mpi = mpi_methods[method_index];

    if (proc_rank == 0)
    {
        my_err = matmul_omp(m1, m2, p_omp, width_omp);
        check(!my_err, "Something went wrong with OpenMP matmul");
    }

    my_err = matmul_mpi(m1, m2, p_mpi, width_mpi,
            proc_rank, num_procs);
    check(!my_err, "Something went wrong with MPI matmul");

    if(m1)
    {
        free(m1);
        debug_mpi(proc_rank, "Freed m1");
    }
    if (m2) 
    {
        free(m2);
        debug_mpi(proc_rank, "Freed m2");
    }

    if (proc_rank == 0)
    {
        for (size_t i = 0; i < width * width; i++)
        {
            double prcnt_err = percent_error(p_omp[i], p_mpi[i]);
            check(prcnt_err < THRESHOLD, "MPI and OMP matmul differs at row %lu, col %lu: %lf %lf", i/width, i%width, p_omp[i], p_mpi[i]);
    
            double elmt;
            scan_rv = scanf("%lf", &elmt);
            check(scan_rv != EOF, "Unexpected EOF");
            check(scan_rv > 0, "Couldn't scan valid double");
    
            prcnt_err = percent_error(p_omp[i], elmt);
            check(prcnt_err < THRESHOLD, "OMP and actual matmul differs at\
                    row %lu, col %lu: %lf %lf", 
                    i/width, i%width, p_omp[i], elmt);
    
            prcnt_err = percent_error(p_mpi[i], elmt);
            check(prcnt_err < THRESHOLD, "MPI and actual matmul differs at\
                    row %lu, col %lu: %lf %lf", 
                    i/width, i%width, p_mpi[i], elmt);
        }
    }
    
    if (proc_rank == 0)
    {
        my_err = gaussian_elimination_naive_inplace_omp(p_omp, width_omp);
        check(!my_err, "Something went wrong during OMP gauss elim");
    }

    my_err = gaussian_elimination_naive_inplace_mpi(p_mpi, width_mpi, 
            proc_rank, num_procs);
    check(!my_err, "Something went wrong during MPI gauss elim");

    if (proc_rank == 0)
    {
        for (size_t i = 0; i < width * width; i++)
        {
            double prcnt_err = percent_error(p_omp[i], p_mpi[i]);            
            debug("%lu %lf %lf", i, p_omp[i], p_mpi[i]);
            check(prcnt_err < THRESHOLD, "Bad gausselim: %lu %lf %lf", i, p_omp[i], p_mpi[i]);
        }
    }

    if (proc_rank == 0)
    {
        log_info("All test cases passed");
    }

    if (p_omp) {
        free(p_omp);
        debug_mpi(proc_rank, "Freed p_omp");
    }
    if (p_mpi) 
    { 
        free(p_mpi);
        debug_mpi(proc_rank, "Freed p_mpi");
    }


    mpi_err = MPI_Finalize();
    debug("MPI environment Finalized");
    check(!mpi_err, "MPI failed to finalize");
    debug_mpi(proc_rank, "Exit success");

    return 0;
error:
    log_warn("Error state reached by process %d", proc_rank);
    if (m1)
        free(m1);
    if (m2)
        free(m2);
    if (p_omp)
        free(p_omp);
    if (p_mpi)
        free(p_mpi);
    mpi_err = MPI_Initialized(&mpi_init_flag);
    if (mpi_err)
        log_warn("Call to `MPI_Initialized` returned with error");
    else if (mpi_init_flag)
    {
        log_info("An MPI environment was found to be initialized.\
                 Finalizing...");
        mpi_err = MPI_Finalize();
        if (mpi_err)
            log_warn("Call to `MPI_Finalize` returned with error");
        else
            log_info("done");
    }
     
}
