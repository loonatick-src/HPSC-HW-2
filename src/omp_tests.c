#include "dbg.h"
#include "impl_omp.h"
#include "matrixio.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>      
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>


impl_omp_t omp_matmul_methods[] = {matMulSquare_baseline_omp,
                                matMulSquare_transpose_omp,
                                matMulSquare_pretranspose_omp};


const int num_methods_omp = sizeof(omp_matmul_methods)/sizeof(impl_omp_t);

const double THRESHOLD = 0.01l;

static
inline
double percent_error(double of, double against)
{
    return fabs(of - against)/against * 100;
}


int main(int argc, char *argv[])
{
    double *m1 = NULL, *m2 = NULL, *p = NULL;
    int widthi, method_index = 1;  // default method_index (transpose)
    
    impl_omp_t matmul = omp_matmul_methods[method_index];
    int num_threads = omp_get_max_threads();

    int scan_rv = scanf("%d", &widthi); 
    check(scan_rv != EOF, "Unexpected EOF");
    check(scan_rv > 0, "Nothing was scanned");

    const int mat_size = widthi * widthi;
    const int mem_size = mat_size * sizeof(double);
    m1 = (double *)malloc(mem_size);
    m2 = (double *)malloc(mem_size);


    int my_err = read_matrices(m1, m2, widthi, stdin);
    check_mem(m1);check_mem(m2);
    check(!my_err, "Something went wrong while reading matrices");

    uint32_t width = (uint32_t) widthi;
    check(width > 0, "There has been an overflow");

    debug("argc: %d", argc);
    if (argc > 1)
    {
        method_index = strtol(argv[1], NULL, 10);
        debug("method index: %d", method_index);
        if (method_index >= num_methods_omp)
        {
            log_warn("Invalid implementation number, proceeding with default");
            method_index = 0;
        }
        matmul = omp_matmul_methods[method_index];
    }

    p = (double *)malloc(width * width * sizeof(double));
    double start_time = omp_get_wtime();
    my_err = matmul(m1, m2, p, width);
    double end_time = omp_get_wtime();
    check(!my_err, "Something went wrong during matrix multiplication");
    double execution_time_matmul = end_time - start_time;
    free(m1);  
    free(m2);

    for (uint32_t i = 0; i < width * width; i++)
    {
        double elmt; 

        scan_rv = scanf("%lf", &elmt);
        check(scan_rv != EOF, "Unexpected EOF");
        check(scan_rv > 0, "Nothing was scanned");
        double prcnt_err = percent_error(p[i], elmt);
        check(prcnt_err < THRESHOLD, "Bad numericals -\
                matrix multiplication test case failed");
    }
    debug("Matrix multiplication complete");
    
    start_time = omp_get_wtime();
    my_err = gaussian_elimination_naive_inplace_omp(p, width);
    end_time = omp_get_wtime();
    check(!my_err, "Something went wrong during gaussian elimination");
    double execution_time_elimination = end_time - start_time;

    /*
    for (int row = 0; row < width; row++)
    {
        for (int col = 0; col < width; col++)
        {
            printf("%lf ", p[row*width + col]);
        }
        putchar('\n');
    }
    */

    debug("Testing complete");
    // width num_threads matmul_time elimination_time
    printf("%u %d %lf %lf\n", width, num_threads,
            execution_time_matmul, execution_time_elimination);
    free(p);
    return 0;
error:
    if (m1) free(m1);
    if (m2) free(m2);
    if (p) free(p);

    return -1;
}
