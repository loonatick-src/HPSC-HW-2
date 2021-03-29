#include "dbg.h"
#include "omp_impl.h"
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


int
main(int argc, char *argv[]) {
    double *M_1 = NULL, *M_2 = NULL, *P = NULL;
    check(argc > 1, "First argument must be width of matrices"); 

    

error:
    if (M_1)
        free(M_1);
    if (M_2)
        free(M_2);
    if (P)
        free(P);
}
