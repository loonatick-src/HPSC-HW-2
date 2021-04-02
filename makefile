CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic
OMPFLAGS=-fopenmp
TESTFLAGS=-pg
SRCDIR=src
BINDIR=bin
TESTDIR=test

omp: $(SRCDIR)/main_omp.c $(SRCDIR)/omp_impl.o
	if [ ! -e $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) -DNDEBUG $(SRCDIR)/main_omp.c $(BINDIR)/omp_impl.o -o $(BINDIR)/omp.out

omp_impl.o: $(SRCDIR)/omp_impl.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -DNDEBUG -c $(SRCDIR)/omp_impl.c -o $(BINDIR)/omp_impl.o

test: test_omp test_mpi
	cd $(TESTDIR) && ./test.out
	cd ..

test_omp: $(SRCDIR)/test.c $(SRCDIR)/omp_impl.c
	if [ ! -e $(TESTDIR) ]; then mkdir $(TESTDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) $(TESTFLAGS) $(SRCDIR)/test.c $(SRCDIR)/omp_impl.c -o $(TESTDIR)/test.out

test_mpi: $(SRCDIR)/test.c
	echo "MPI testing not implemented"

clean:
	rm $(BINDIR)/*.out
	rm $(BINDIR)/*.o
	rm $(TESTDIR)/*.out
	rmdir $(BINDIR)
