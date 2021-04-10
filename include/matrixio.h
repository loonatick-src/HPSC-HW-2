#ifndef _MATRIX_IO_H
#define _MATRIX_IO_H
/* Header only */

#include <stdio.h>
#include <stdlib.h>
int
read_matrix_intonull(FILE *file, double *m, int width)
{
    // assumes m is NULL
    check(m == NULL, "Expects a NULL pointer, initialized pointer passed");
    check(width > 0, "Width < 0");
    const int mat_size = width * width;
    const int mem_size = mat_size * sizeof(double);
    m = (double *)malloc(mem_size);
    check_mem(m);

    for (int i = 0; i < mat_size; i++)
    {
        int scan_flag = fscanf(file, "%lf", m+i);
        check(scan_flag != EOF, "Unexpected EOF");
        check(scan_flag > 0, "Nothing was scanned");
    }

    return 0;
error:
    return 10;
}


int
read_matrix(FILE* file, double *m, int width)
{
    check_mem(m);
    check(width > 0, "Width < 0");
    const int mat_size = width * width;
    for (int i = 0 ; i < mat_size; i++)
    {
        fscanf(file, "%lf", m+i);
        check(scan_flag != EOF, "Unexpected EOF");
        check(scan_flag > 0, "Nothing was scanned");
    }
    return 0;
error:
    // let caller handle error
    return -1;
}


int
read_matrices(double *m1, double *m2, int *width_p, FILE *file)
{
    check(m1 == NULL && m2 == NULL, "NULL pointers expected");
    check(file, "Not a valid FILE pointer");
    int scan_flag = fscanf(file, "%d", width_p);
    check(scan_flag != EOF, "Unexpected EOF");
    check(scan_flag > 0, "Nothing was scanned");
    check(*width > 0, "`width` appears to be a negative number.");

    int my_err = read_matrix_intonull(file, m1, width);
    check(my_err, "Something went wrong while reading first matrix");
    my_err = read_matrix_intonull(file, m2, width);
    check(my_err, "Something went wrong while reading second matrix");

    return EXIT_SUCCESS;
error:
    // let caller handle error
    return EXIT_FAILURE;
}


int
print_matrix(FILE *file, double *m, int width)
{
    check(file != NULL, "file is NULL");
    check_mem(m);
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < width-1; j++)
        {
            print_flag = fprintf(file, "%lf ", m[i*width + j]);
        }
        fprintf(file, "%lf\n", m[i*width + width - 1]);
    }
    return EXIT_SUCCESS;
error:
    return EXIT_FAILURE;
}

#endif
