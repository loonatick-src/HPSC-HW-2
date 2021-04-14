#include "dbg.h"
#include "impl_omp.h"

#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>

#define mel(A, w, i, j) A[i*w + j]  // get matrix element row,col
#define submel(A, w, s, r, c) A[(r+s)*w + c + s]  // get principal submatrix element at row,col

#define SGN64(A) ((0x9000000000000000&(A))>>63)

typedef struct {
    double g;
    double s;
} givens_t;


int
matMulSquare_baseline_omp(const double *M_1,
                      const double *M_2,
                      double *P, 
                      const uint32_t width) 
{
    check_mem(M_1); check_mem(M_2); check_mem(P);
    debug("Performing matrix multiplication for %d x %d matrices", width, width);
    const uint32_t matrix_size = width * width;
    check(matrix_size >= width, "Integer overflow (uint32_t).");
    debug("%u %u", width, matrix_size);
#   pragma omp parallel for
    for (uint32_t row = 0; row < width; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            const uint32_t i_p = row * width + col;
            double elmt_sum = 0.0l;
            // OpenMP implementation parallelizes the calculation of
            // individual matrix elements of the product matrix
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
matMulSquare_transpose_omp(const double *M_1,
                       const double *M_2,
                       double *P,
                       uint32_t width)
{
    // transposing the second matrix
    // fewer cache misses
    // large overhead of transposition
    // cumulative less than baseline?
    double *M_2trnsps = NULL;
    const uint32_t matrix_size = width * width;
    check(matrix_size >= width, "Integer overflow (uint32_t).");

    const unsigned long int mem_size = matrix_size * sizeof(double);
    check(mem_size >= matrix_size, "Integer overflow when calculating size to be malloc'd");

    // 1m
    M_2trnsps = (double *) malloc(mem_size);
    check_mem(M_2trnsps);

    // parallelized matrix transposition
#   pragma omp parallel for
    for (uint32_t row = 0; row < width; row++ )
    {
        for (uint32_t col = 0; col < width; col++)
        {
            M_2trnsps[col * width + row] = M_2[row * width + col]; 
        }
    }

    // matrix multiplication, taking into account transposition
    // of M_2
#   pragma omp parallel for
    for (uint32_t row = 0; row < width; row++) 
    {
        for (uint32_t col = 0; col < width; col++)
        {
            const uint32_t i_p = row * width + col;
            double elmt_sum = 0.0l;
            for (uint32_t i = 0; i < width; i++)
            {
                elmt_sum += M_1[row*width + i]*M_2trnsps[col * width + i];
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
matMulSquare_pretranspose_omp(const double *M_1,
                          const double *M_2,
                          double *P,
                          uint32_t width)
{
    const uint32_t matrix_size = width * width;
    check(matrix_size >= width, "Integer overflow (uint32_t).");
    
#   pragma omp parallel for
    for (uint32_t row = 0; row < width; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            const uint32_t i_p = row * width + col;
            double elmt_sum = 0.0l;
            for (uint32_t i = 0; i < width; i++)
            {
                elmt_sum += M_1[row*width + i] * M_2[col*width + i];
            }
            P[i_p] = elmt_sum;
        }
    }

    return 0;
error:
    return -1;
}


int
gaussian_elimination_naive_inplace_omp(double *M, uint32_t width)
{
    check_mem(M);
    for (uint32_t iter = 0; iter < width - 1; iter++) 
    {
        double pivot = M[iter * width + iter];
        check(pivot != 0, "Zero pivot found! Use partial pivoting algo.");
#       pragma omp parallel for
        for (uint32_t row = iter+1; row < width; row++)
        {
            double leverage = M[row*width + iter];  // LOL
            for (uint32_t col = iter + 1; col < width; col++)
            {
                M[row*width + col] -= leverage/pivot;
            }
            M[row * width + iter] = 0.0l;
        }
    }
    return 0;
error:
    return -1;
}
