#ifndef _MATRIX_IO_H
#define _MATRIX_IO_H
/* Header only */

#include <stdio.h>
#include <stdlib.h>

int
read_matrix(FILE* file, double *m, int width)
{
    check_mem(m);
    check(width > 0, "Width < 0");
    const int mat_size = width * width;
    for (int i = 0 ; i < mat_size; i++)
    {
        int scan_flag = fscanf(file, "%lf", m+i);
        check(scan_flag != EOF, "Unexpected EOF");
        check(scan_flag > 0, "Nothing was scanned");
    }
    return 0;
error:
    // let caller handle error
    return -1;
}


int
read_matrices(double *m1, double *m2, int width, FILE *file)
{
    check_mem(m1);check_mem(m2);
    check(file, "Not a valid FILE pointer");
    check(width > 0, "Non-positive width");

    int my_err = read_matrix(file, m1, width);
    check_mem(m1);
    my_err = read_matrix(file, m2, width);
    check(!my_err, "Something went wrong while reading second matrix");
    check_mem(m2);

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
            fprintf(file, "%lf ", m[i*width + j]);
        }
        fprintf(file, "%lf\n", m[i*width + width - 1]);
    }
    return EXIT_SUCCESS;
error:
    return EXIT_FAILURE;
}

#endif
