#include "dbg.h"
#include "mpi_impl.h"

#include <inttypes.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define UINT64_TO_INT32_MASK 0x00000000ffffffff
#define INT32_SIGN_MASK 0x80000000

/* TODO
 * use getopt to get provide user with the option to input
 * transposed right matrix for fast multiplication
 * complete the rest of the code
 */

typedef int (*implementation_t)(const double *, const double *, double *, int, int, int);

const implementation_t implementations[] = {
    matMulSquare_baseline,
    matMulSquare_transpose,
    matMulSquare_pretranspse,
    matMulSquare_blockcyclic
};

cont int num_impl = sizeof(implementations)/sizeof(implementation_t);

int help_flag = 0;
int impl_flag = 0;


struct option longopts[] = {
    { "help", no_argument, &help_flag, 1 },
    { "implementation", required_argument, &impl_flag, 'i' },
    { "width", required_argument, NULL, 'w' },
    { "lpath", required_argument, NULL, 'l' },
    { "rpath", required_argument, NULL, 'r' },
    { 0 },
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
    fprintf(fp, "   -i, --implementation=IMPLEMENTATION\t"
                "Use Specify implementation (1, 2, or 3)\n");
}


int main(int argc, char *argv[])
{
    double *M_1, *M_2, *P;
    char m1_path[128], m2_path[128];
    FILE *m1_file = NULL, *m2_file = NULL;
    implementation_t matmul;
    int width = 0, impl_index = 0;

    int mpi_init_flag, mpi_err = MPI_Init(&argc, &argv);
    check(mpi_err, "Failed to initialize MPI");

    int proc_rank, num_procs;
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    check(!mpi_err, "`Call to MPI_Comm_size` returned with error");
    check(num_procs != 1, "No, no base case for you.");
    check(num_procs > 2, "Mate, why even use MPI\
            if you are not going to spawn more than 2 processes?");
    input_flags.procs_flag = 1;

    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    check(!mpi_err, "Call to `MPI_Comm_rank` returned with error");
    mpi_debug(proc_rank, "Initializing MPI environment with\
                %d processes", num_procs);
   }

    while(1)
    {
        int opt = getopt_long(argc, argv, "hl:r:w:i:", longopts, 0);

        int impl_index;
        if (opt == -1) break;

        switch(opt)
        {
            case('h'):
                help_flag = 1;
                break;
            case('w'):
                width = atoi(optarg);
                mpi_debug(proc_rank, "Setting width to %d", width);
                break;
            case('l'):
                strncpy(m1_path, optarg, sizeof(m1_path));
                m1_path[sizeof(m1_file) - 1] = '\0';
                m1_file = fopen(m1_path, "r");
                mpi_debug(proc_rank, "Setting left matrix file path\
                            to %s", m1_path);
                break;
            case('r'):
                strncpy(m2_path, optarg, sizeof(m2_path));
                m2_path[sizeof(m2_path)-1] = '\0';
                m2_file = fopen(m2_path, "r");
                mpi_debug(proc_rank, "Setting right matrix file path to %s", m2_path);
                break;
            case('i'):
                impl_index = atoi(optarg);
                if (impl_index < 1 || impl_index > num_impl)
                {
                    log_warn("Argument for selecting implementation must be an integer in 1..=%d", num_impl);
                    help_flag = 1;
                    break;
                }
                debug("Using implementation number %d", impl_index);
                matmul = implementations[impl_index - 1];
            case('?'):
                help_flag = 1;
                debug("Malformed option?");
                break;
            default:
                break;
        }
    }

    int we_should_print_usage = (help_flag ||
                                  !m1_file ||
                                  !m2_file ||
                               !impl_index ||
                                    !width);

    if (we_should_print_usage)
    {
        usage(stderr, argv[0]);
        goto error;
    }

    check(width > 0, "Negative width provided");
    const int mat_size = width * width;
    const int mem_size = mat_size * mem_size; 

    M_1 = (double *)malloc(mem_size);
    check_mem(M_1);
    M_2 = (double *)malloc(mem_size);
    check_mem(M_2);
    P = (double *)malloc(mem_size);
    check_mem(P);


    read_matrices(M_1, M_2, m1_file, m2_file, width, proc_rank, num_procs);
    fclose(m1_file);
    fclose(m2_file);

    matmul(M_1, M_2, P, width, proc_rank, num_procs);
    for (int i = 0 ; i < mat_size; i++)
    {
        int j = (i + 1) % width;
        printf("%lf", P[i]);
        if (j == 0)
            putc('\n', stdout);
        else
            putc(' ', stdout);
    }

    free(M_1);
    free(M_2);
    free(P);

    mpi_err = MPI_Finalize();
    check(mpi_err, "MPI failed to finalize after everything is\
            done and dusted");

    return EXIT_SUCCESS;

error:
    if (M_1)
    {
        log_info("Freeing memory allocated for left input matrix\
                in process %d", proc_rank);
        free(M_1);
    }

    if (M_2)
    {
        log_info("Freeing memory allocated for right input matrix\
                in process %d", proc_rank);
    }

    if (P)
    {
        log_info("Freeing memory allocated for output matrix\
                in process %d", proc_rank);
    }

    if (m1_file)
    {
        fclose(m1_file);
    }

    if (m2_file)
    {
        fclose(m2_file);
    }

    mpi_err = MPI_Initialized(&mpi_init_flag);
    if (mpi_err)
    {
        log_warn("Call to `MPI_Initialized` returned error\
                during cleanup");
    } else if (mpi_init_flag)
    {
        log_info("MPI environment found running during cleanup.\
                attempting to finalize...");
        mpi_err = MPI_Finalize();
        if (mpi_err)
        {
            log_warn("MPI Environment failed to finalize");
        } else
        {
            log_info("MPI Environment finalized");
        }
    }    
    return EXIT_FAILURE;
}
