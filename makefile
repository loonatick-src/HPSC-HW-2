CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic
OMPFLAGS=-fopenmp
TESTFLAGS=-pg
SRCDIR=src
BINDIR=bin
TESTDIR=test

omp: main_omp.c omp_impl.o
	if [ ! -e $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) -DNDEBUG $(SRCDIR)/main_omp.c $(BINDIR)/omp_impl.o -o $(BINDIR)/omp.out

omp_impl.o: omp_impl.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -DNDEBUG -c $(SRCDIR)/omp_impl.c -o $(BINDIR)/omp_impl.o

test: test_omp test_mpi

test_omp: main_omp.c omp_impl.c
	if [ ! -e $(TESTDIR) ]; then mkdir $(TESTDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) $(TESTFLAGS) $(SRCDIR)/omp_impl.c $(SRCDIR)/main_omp.c -o $(TESTDIR)/omp_test.out

test_mpi: main_mpi.c
	echo "MPI testing not implemented"

clean:
	rm ../bin/*
	rmdir ../bin
