#include "dbg.h"
#include "omp_impl.h"

#include <getopt.h>
#include <inttypes.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef int (*implementation_t)(const double *, const double *, double *, uint32_t);


typedef struct {
    char lmat_flag;
    char rmat_flag;
    char w_flag;
    char thread_flag;
} flags_t;


flags_t input_flags = { 0, 0, 0, 0 };


const implementation_t implementations[] = {
    matMulSquare_baseline,
    matMulSquare_transpose,
    matMulSquare_pretranspose };
const int num_impl = sizeof(implementations) / sizeof(implementation_t);

int help_flag = 0;
int impl_flag = 0;


struct option longopts[] = {
    { "help", no_argument, &help_flag, 1 },
    { "implementation", required_argument, &impl_flag, 'i' },
    { "width", required_argument, NULL, 'w'},
    { "lpath", required_argument, NULL, 'l'},
    { "rpath", required_argument, NULL, 'r'},
    { "threadcount", required_argument, NULL, 't'},
    { 0 }
};



void
usage(FILE *fp, const char *path)
{
    // TODO: complete this help message
    const char *basename = strchr(path, '/');
    basename = basename ? basename + 1 : path;

    fprintf(fp, "usage: %s [OPTION]\n", basename);
    fprintf(fp, "   -h, --help\t\t"
                "Print this help and exit.\n");
    fprintf(fp, "   -w, --width=WIDTH\t"
                "Width of square matrices to be input\n");
    fprintf(fp, "   -l, --lpath=LEFT_MATRIX_FILE_PATH\t"
                "Path to (row-major) entries of left matrix.\n");
    fprintf(fp, "   -r, --rightpath=RIGHT_MATRIX_FILE_PATH\t"
                "Path to (row-major) entries of right matrix.\n");
    fprintf(fp, "   -t, --threadcount=OMP_NUM_THREADS\t"
                "Number of threads to be used by OpenMP.\n");
    fprintf(fp, "   -i, --implementation=IMPLEMENTATION\n"
                "Use Specify implementation (1, 2, or 3)\n");
}



int
main(int argc, char *argv[])
{
    double *M_1 = NULL, *M_2 = NULL, *P = NULL;
    char m1_path[128] = { 0 }, m2_path[128] = { 0 };  // fixed size buffers
    FILE *m1_file = NULL, *m2_file = NULL;
    int widthi;
    implementation_t mulfunc; 

    if (argc < 9)
    {
        usage(stderr, argv[0]);
        goto error;
    }
    debug("Processing inputs");
    while(1)
    {
        int opt = getopt_long(argc, argv, "hl:r:t:w:i:", longopts, 0);
        int index;
        if (opt == -1) break;

        switch(opt)
        {
            case('h'):
                help_flag = 1;
                break;
            case('w'):
                widthi = atoi(optarg);
                input_flags.w_flag = 1;
                debug("Setting width to %d", widthi);
                break;
            case('i'):
                index = atoi(optarg);
                if (index < 1 || index > num_impl)
                {
                    log_err("Argument for selecting implementations must be an integer in 1..=%d", num_impl);
                    help_flag = 1;
                    break;
                }
                mulfunc = implementations[index-1];
                debug("Using implementation number %d", index);
                impl_flag = 1;
                break;
            case('l'):
                input_flags.lmat_flag = 1;
                strncpy(m1_path, optarg, sizeof(m1_file));
                m1_path[sizeof(m1_file) - 1] = '\0';
                debug("Left matrix filename: %s", m1_path);
                break;
            case('r'):
                strncpy(m2_path, optarg, sizeof(m2_file));
                m2_path[sizeof(m2_path) - 1] = '\0';
                input_flags.rmat_flag = 1;
                debug("Right matrix file name: %s", m2_path);
                break;
            case('t'):
                // no return value, I hope that its error handling is robust
                input_flags.thread_flag = 1;
                debug("Setting number of threads to %d", atoi(optarg));
                omp_set_num_threads(atoi(optarg));
                break;
            case('?'):
                help_flag = 1;
                debug("Malformed argument?");
                break;
            default:
                break;
        }
    }


    int we_should_print_usage = (help_flag ||
                    !input_flags.rmat_flag ||
                    !input_flags.lmat_flag ||
                  !input_flags.thread_flag ||
                                !impl_flag);
    if (we_should_print_usage)
    {
        usage(stderr, argv[0]);
        goto error;
    }
    debug("input processed successfully");


    // all matrices width x width
    check(widthi > 0, "Width of matrices must be a postive integer\
            (provided value: %d", widthi);
    const uint32_t width = (uint32_t) widthi;
    debug("Setting width to %u", width);

    // number of elements per matrix
    const uint32_t matsize = width * width;
    // check for integer overflow
    check(matsize >= width, "Proposed matrix size too large\
            (integer overflow)");

    // allocating memory for storing matrices
    const uint32_t mem_size = matsize * sizeof(double);
    // check for integer overflow
    check(mem_size > matsize, "Proposed matrix size too large\
            (integer overflow while allocating memory)");

    // malloc
    // 1m
    debug("Allocating memory to M1");
    M_1 = (double *) malloc(mem_size);
    check_mem(M_1);

    // 2m
    debug("Allocating memory to M2");
    M_2 = (double *) malloc(mem_size);
    check_mem(M_2);

    // 3m
    debug("Allocating memory to P");
    P = (double *) malloc(mem_size);
    check_mem(P);


    // reading matrix elements from file
    m1_file = fopen(m1_path, "r");
    debug("Opening file for reading M1");

    uint32_t in_count = 0;
    double melmnt;
    while (in_count < matsize)
    {
        int scanf_err = fscanf(m1_file, "%lf", &melmnt);
        check(scanf_err != EOF, "Unexpected EOF,\
                insufficient  matrix elements provided\
                (Promised %u, provided %u", matsize, in_count);
        check(scanf_err, "`scanf` returned with error\
                after reading %u matrix elements", in_count);
        M_1[in_count] = melmnt;
        in_count++;
    }
    fclose(m1_file);
    debug("M1 read complete. Moving on to M2");

    m2_file = fopen(m2_path, "r");
    in_count = 0;
    while (in_count < matsize)
    {
        int scanf_err = fscanf(m2_file, "%lf", &melmnt);
        check(scanf_err != EOF, "Unexpected EOF,\
                insufficient  matrix elements provided\
                (Promised %u, provided %u", matsize, in_count);
        check(scanf_err, "`scanf` returned with error\
                after reading %u matrix elements", in_count);
        M_2[in_count] = melmnt;
        in_count++;
    }
    fclose(m2_file);
    debug("M2 read complete");


    debug("Performing matrix multiplication");
    int my_err = mulfunc(M_1, M_2, P, width);
    check(my_err, "Something went wrong during matrix multiplication");

    debug("Returned from matrix multiplication");
    debug("Printing product matrix to stdout");
    for (int row = 0; row < width; row++)
    {
        int skip = width * row;
        for (int col = 0; col < width-1; col++)
        {
            printf("%lf ", P[skip + col]);
        }
        printf("%lf\n", P[skip + width-1]);
    }

     
    free(M_1);  // 1m
    free(M_2);  // 2m
    free(P);    // 3m

    return EXIT_SUCCESS;

error:
    if (M_1)
        free(M_1);  // 1m
    if (M_2)
        free(M_2);  // 2m
    if (P)
        free(P);    // 3m
    if (m1_file)
        fclose(m1_file);
    if (m2_file)
        fclose(m2_file);
    return EXIT_FAILURE;
}
