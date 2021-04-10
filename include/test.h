#ifndef __TEST__H_
#define __TEST__H_
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#define threshold_percent_error 1.0e-6

double
percent_error(double of, double against)
{
    return (fabs(of - against)/against) * 100.0l;
}


int
validate_matrix(double *v, double *a, uint32_t matsize)
{
    check_mem(v);
    check_mem(a);
    for (uint32_t index = 0; index < matsize; index++)
    {
        double err = percent_error(v[index], a[index]);
        if (err > threshold_percent_error)
        {
            log_err("Erroneous value at index %d, expected close to %lf, found %lf", index, a[index], v[index]);
            return -1;
        }
        index++;
    }
    return 0;
}
#endif
