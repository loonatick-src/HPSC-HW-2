#include "dbg.h"
#include "omp_impl.h"
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>


void
matMulSquare_baseline(const double * const M_1,
                      const double * const M_2,
                      double *P, 
                      const uint32_t width) 
{
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
}


int
matMulSquare_transpose(const double * const M_1,
                       const double * const M_2,
                       double *P,
                       const uint32_t width)
{
    double *M_2trnsps = NULL;
    const uint32_t matrix_size = width * width;
    check(matrix_size >= width, "Integer overflow (uint32_t).");

    const unsigned long int mem_size = matrix_size * sizeof(uint32_t);
    check(mem_size >= matrix_size, "Integer overflow when calculating size to be malloc'd");

    M_2trnsps = (double *) malloc(mem_size);
    mem_check(M_2trnsps);

    for (uint32_t row = 0; row < width; row++ )
    {
#       pragma omp parallel for
        for (uint32_t col = 0; col < width; col++)
        {
            M_2trnsps[col * width + row] = M_2[row * width + col]; 
        }
    }

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
                sum += M_1[index] * M_2trnsps[index];
            }
            P[i_p] = sum;
        }
    }

    free(M_2trnsps);
    return 0;

error:
    if (M_2trnsps)
        free(M_2trnsps);
    return -1;
}
