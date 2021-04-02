#ifndef _MY_DEBUG_H
#define _MY_DEBUG_H

/* Useful debugging macros
 * taken directly from Zed Shaw's
 * Learn C the Hard Way */

#include <stdio.h>
#include <errno.h>
#include <string.h>

#ifdef NDEBUG
#define debug(M, ...)
#define mpi_debug(R, M, ...)
#else
#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define mpi_debug(R, M, ...) if ((R) == 0) { debug(M, ##__VA_ARGS__); }
#endif

#define clean_errno() (errno == 0 ? "None" : strerror(errno))

#define log_err(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

#define log_warn(M, ...) fprintf(stderr, "[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

#define log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n", __FILE__, __LINE__ ##__VA_ARGS__)

#define check(A, M, ...) if(!(A)) { log_err(M, ##__VA_ARGS__); errno=0; goto error;}

#define sentinel(M, ...) { log_err(M, ##__VA_ARGS__); errno=0; gotoerror; }

#define check_mem(A) check((A), "Out of memory.")

#define check_debug(A, M, ...) if (!(A)) { debug(M, ##__VA_ARGS__); errno=0; goto error; }


#endif
