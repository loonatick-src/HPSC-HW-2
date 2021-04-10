CFLAGS = -Wall -Winline -pedantic -pg -Iinclude
LDLIBS = -lm

mpi:
	mpicc $(CFLAGS) -DBALANCED_BASE src/impl_mpi.c src/test_mpi.c -o bin/mpi.out $(LDLIBS)
