mpi: CFLAGS = -Wall -Wextra -pedantic -O2 -Iinclude -DNDEBUG
mpi:
	mpicc $(CFLAGS) src/impl_mpi.c src/mpi_tests.c -lm -o bin/mpi.out 

omp: CFLAGS = -fopenmp -Wall -Wextra -pedantic -O2 -Iinclude -DNDEBUG
omp:
	gcc $(CFLAGS) src/impl_omp.c src/omp_tests.c -lm -o bin/omp.out

