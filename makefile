mpi: CFLAGS = -Wall -Wextra -Werror -pedantic -O2 -Iinclude
mpi:
	mpicc $(CFLAGS) src/impl_mpi.c src/mpi_tests.c -lm -o bin/mpi.out 

omp: CFLAGS = -fopenmp -Wall -Wextra -Werror -pedantic -O2 -Iinclude
omp:
	gcc $(CFLAGS) src/impl_omp.c src/omp_tests.c -lm -o bin/omp.out

gelim: CFLAGS = -fopenmp -Wall -Wextra -Werror -pedantic -O2 -Iinclude
gelim:
	mpicc $(CFLAGS) src/impl_omp.c src/impl_mpi.c src/test_gelim.c -lm -o bin/gelim.out
