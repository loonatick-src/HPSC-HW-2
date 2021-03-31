#include "dbg.h"
#include "omp_impl.h"

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
const double threshold_percent_error = 1.e-4;
const char *m1_path = "m_1_100.dat";
const char *m2_path = "m_2_100.dat";
const char *m2t_path = "m_2t_100.dat";
const char *mmulpath = "matmul_100.dat";
const char *mmulpatht = "matmul_100t.dat";


int
read_matrix(double *m, FILE *filehandle, uint32_t matsize)
{
    check(matsize > 0, "Matrix size must be a positive integer, duh");
    for (int i = 0; i < matsize; i++)
    {
        int fscanf_err = fscanf(filehandle, "%lf", &(m[i]));
        check(fscanf_err, "fscanf threw some error");
    }
    return EXIT_SUCCESS;
error:
    return EXIT_FAILURE;
    
}


double
percent_error(double of, double against)
{
    return (fabs(of - against)/against) * 100.0l;
}


int
validate_matrix(double *v, double *a, uint32_t matsize)
{
    for (uint32_t index = 0; index < matsize; index++)
    {
        double err = percent_error(v[index], a[index]);
        if (err > threshold_percent_error)
            return -1;
        index++;
    }
    return 0;
}


void
test_matmul_omp()
{
    log_info("Testing OpenMP matrix multiplication implementations");
    log_info("Allocating space for matrices");
    double *m1 = (double *)malloc(memsize);
    double *m2 = (double *)malloc(memsize);
    double *matmul = (double *)malloc(memsize);
    double *p = (double *)malloc(memsize);

    FILE *m1_file = NULL;
    FILE *m2_file = NULL;
    FILE *m2t_file = NULL;
    FILE *matmul_file = NULL;

    log_info("Opening files for reading matrices");
    m1_file = fopen(m1_path, "r");
    m2_file = fopen(m2_path, "r");
    m2t_file = fopen(m2t_path, "r");
    matmul_file = fopen(mmulpath, "r");

    log_info("Reading matrices into memory");
    read_matrix(m1, m1_file, matsize);
    fclose(m1_file);

    read_matrix(m2, m2_file, matsize);
    fclose(m2_file);

    read_matrix(matmul, matmul_file, matsize);
    fclose(matmul_file);


    log_info("Testing first implementation (direct multiplication)");
    int my_err = matMulSquare_baseline(m1, m2, p, width);
    check(my_err, "Something went wrong in the first implementation");
    my_err = validate_matrix(p, matmul, matsize);
    check(my_err == 0, "Output of first implementation seems erroneous");

    log_info("Testing second implementation (transposed multiplication)");
    my_err = matMulSquare_transpose(m1, m2, p, width);
    check(!my_err, "Something went wrong in the second implementation");


    log_info("Reading transposed matrix into memory");
    my_err = read_matrix(m2, m2t_file, matsize);
    check(my_err, "Something went wrong while reading transposed matrix");
    fclose(m2t_path);

    log_info("Testing second implementation (pretransposed multiplication");
    my_err = matMulSquare_pretranspose(m1, m2, p, width);

    free(m1);
    free(m2);
    free(p);
    free(matmul);
error:
    if (m1) free(m1);
    if (m2) free(m2);
    if (p) free(p);


}
