CFLAGS = -Wall -Winline -pedantic -g -Iinclude
LDLIBS = -lm

mpi: CFLAGS = -Wall -Winline -pedantic -g -Iinclude
mpi:
	mpicc $(CFLAGS) src/impl_mpi.c src/test_mpi.c -o bin/mpi.out $(LDLIBS)

omp: CFLAGS = -Wall -Winline -pedantic -g -Iinclude -fopenmp
omp:
	gcc $(CFLAGS) src/impl_omp.c src/test_omp.c -o bin/omp.out $(LDLIBS)

