#include "dbg.h"
#include "omp_impl.h"

#include <inttypes.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>


int
matMulSquare_baseline(const double *M_1,
                      const double *M_2,
                      double *P, 
                      const uint32_t width) 
{
    debug("Performing matrix multiplication for %d x %d matrices", width, width);
    const uint32_t matrix_size = width * width;
    check(matrix_size >= width, "Integer overflow (uint32_t).");
    debug("%u %u", width, matrix_size);
    for (uint32_t row = 0; row < width; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            const uint32_t i_p = row * width + col;
            double elmt_sum = 0.0l;
            // OpenMP implementation parallelizes the calculation of
            // individual matrix elements of the product matrix
#           pragma omp parallel for reduction(+:elmt_sum)
            for (uint32_t i = 0; i < width; i++) 
            {
                elmt_sum += M_1[row*width + i] * M_2[i*width + col];
            }
            P[i_p] = elmt_sum;
        }
    }
    debug("Matrix multiplication complete");

    return 0;
error:
    return -1;
}


int
matMulSquare_transpose(const double *M_1,
                       const double *M_2,
                       double *P,
                       uint32_t width)
{
    // transposing the second matrix
    // fewer cache misses
    // large overhead of transposition
    double *M_2trnsps = NULL;
    const uint32_t matrix_size = width * width;
    check(matrix_size >= width, "Integer overflow (uint32_t).");

    const unsigned long int mem_size = matrix_size * sizeof(double);
    check(mem_size >= matrix_size, "Integer overflow when calculating size to be malloc'd");

    // 1m
    M_2trnsps = (double *) malloc(mem_size);
    check_mem(M_2trnsps);

    // parallelized matrix transposition
    for (uint32_t row = 0; row < width; row++ )
    {
#       pragma omp parallel for
        for (uint32_t col = 0; col < width; col++)
        {
            M_2trnsps[col * width + row] = M_2[row * width + col]; 
        }
    }

    // matrix multiplication, taking into account transposition
    // of M_2
    for (uint32_t row = 0; row < width; row++) 
    {
        for (uint32_t col = 0; col < width; col++)
        {
            const uint32_t i_p = row * width + col;
            double elmt_sum = 0.0l;
#           pragma omp parallel for reduction(+:elmt_sum)
            for (uint32_t i = 0; i < width; i++)
            {
                const uint32_t index = row * width + i;
                elmt_sum += M_1[index] * M_2trnsps[index];
            }
            P[i_p] = elmt_sum;
        }
    }

    free(M_2trnsps);   // 1m
    return 0;

error:
    if (M_2trnsps)
        free(M_2trnsps);
    return -1;
}


int
matMulSquare_pretranspose(const double *M_1,
                          const double *M_2,
                          double *P,
                          uint32_t width)
{
    const uint32_t matrix_size = width * width;
    check(matrix_size >= width, "Integer overflow (uint32_t).");
    
    for (uint32_t row = 0; row < width; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            const uint32_t i_p = row * width + col;
            double elmt_sum = 0.0l;
#           pragma omp parallel for reduction(+:elmt_sum)
            for (uint32_t i = 0; i < width; i++)
            {
                const uint32_t index = row * width + i;
                elmt_sum += M_1[index] * M_2[index];
            }
            P[i_p] = elmt_sum;
        }
    }

    return 0;
error:
    return -1;
}
