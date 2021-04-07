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


int
matMulSquare_baseline_omp(const double *M_1,
                      const double *M_2,
                      double *P, 
                      const uint32_t width) 
{
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


double
norm_l2(const double *vector, uint32_t dim)
{
    double rv = 0.0l;
#   pragma omp parallel for reduction(+:rv)
    for(uint32_t i = 0; i < dim; i++)
    {
        rv += vector[i] * vector[i];
    }
    return (sqrt(rv));
}


double
inner_product(double *v_1, double *v_2, uint32_t dim)
{
    double rv = 0.0l;
#   pragma omp parallel for reduction(+:rv)
    for (uint32_t i = 0; i < dim; i++)
    {
        rv += v_1[i] * v_2[i];
    }

    return rv;
}


int
eye_calloc(double *matrix, uint32_t width)
{
    check_mem(matrix);
#   pragma omp parallel for shared(matrix)
    for (uint32_t i = 0; i < width; i++)
    {
        matrix[i*width + i] = 1.0l;
    }
    return EXIT_SUCCESS;
error:
    return EXIT_FAILURE;
}


int
unit_vector(const double *vector_in, double *vector_out, uint32_t width)
{
    check_mem(vector_in);
    check_mem(vector_out);
    double norm = norm_l2(vector_in, width);
    check(norm > 0, "Null vector passed");
#   pragma omp parallel for shared(vector_out)
    for (uint32_t i = 0; i < width; i++)
    {
        vector_out[i] = vector_in[i]/norm;
    }
    return EXIT_SUCCESS;
error:
    return EXIT_FAILURE;
}


int
unit_vector_inplace(double *vector, uint32_t width)
{
    check_mem(vector);
    double norm = norm_l2(vector, width);
    check(norm > 0, "Null vector passed");
    for (uint32_t i = 0; i < width; i++)
    {
        vector[i]/=norm;
    }
    return EXIT_SUCCESS;
error:
    return EXIT_FAILURE;
}


int
eye_col_calloc(double *vector, uint32_t width, uint32_t n)
{
    check_mem(vector);
    vector[n] = 1.0l;
    return EXIT_SUCCESS;

error:
    return EXIT_FAILURE;
}
