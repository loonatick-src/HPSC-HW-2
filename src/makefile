CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic
OMPFLAGS=-fopenmp
DBGFLAGS=-g
BINDIR=../bin
TESTDIR=../test

omp: main_omp.c omp_impl.o
	if [ ! -e $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) -DNDEBUG main_omp.c omp_impl.o -o $(BINDIR)/omp.out

omp_impl.o: omp_impl.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c omp_impl.c

test: test_omp test_mpi

test_omp: main_omp.c omp_impl.c
	if [ ! -e $(TESTDIR) ]; then mkdir $(TESTDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) $(DBGFLAGS) omp_impl.c main_omp.c -o ../test/omp_test.out

test_mpi: main_mpi.c
	echo "MPI testing not implemented"

clean:
	rm ../bin/*
	rmdir ../bin
