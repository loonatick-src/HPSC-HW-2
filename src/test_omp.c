#include "dbg.h"
#include "impl_omp.h"
#include "matrixio.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>      
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>

implementation_t omp_matmul_methods[] = {matMulSquare_baseline_omp,
                                matMulSquare_transpose_omp,
                                matMulSquare_pretranspose_omp};


const int num_methods_omp = sizeof(omp_matmul_methods)/sizeof(implementation_t);


int main(int argc, char *argv[])
{
    double *m1 = NULL, *m2 = NULL, *p = NULL;
    int widthi, method_index = 1;  // default method_index (transpose)
    
    implementation_t matmul = omp_matmul_methods[method_index];
    omp_set_num_threads(6);

    scanf("%d", &widthi); 
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
    my_err = matmul(m1, m2, p, width);
    check(!my_err, "Something went wrong during matrix multiplication");

    free(m1);  
    free(m2);

    for (int i = 0; i < width * width; i++)
    {
        double elmt; scanf("%lf", &elmt);
        if (fabs(elmt - p[i])/elmt > 0.01)
        {
            log_warn("Bad numericals: expected %lf, found %lf", elmt, p[i]);
        }
    }
    debug("Matrix multiplication complete");
    
    /*
    my_err = gaussian_elimination_naive_inplace(p, width);
    check(!my_err, "Something went wrong during gaussian elimination");

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
    free(p);
    return 0;
error:
    if (m1) free(m1);
    if (m2) free(m2);
    if (p) free(p);

    return -1;
}
