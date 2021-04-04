CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic
OMPFLAGS=-fopenmp
TESTFLAGS=-pg $(CFLAGS)
SRCDIR=src
BINDIR=bin
TESTDIR=test
TESTFILES=("$(TESTDIR)/m_1_100.dat","$(TESTDIR)/m_2_100.dat","$(TESTDIR)/matmul_100.dat","$(TESTDIR)/matmul_100t.dat")

omp: $(SRCDIR)/main_omp.c $(SRCDIR)/omp_impl.o
	if [ ! -e $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) -DNDEBUG $(SRCDIR)/main_omp.c $(BINDIR)/omp_impl.o -o $(BINDIR)/omp.out


omp_impl.o: $(SRCDIR)/omp_impl.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -DNDEBUG -c $(SRCDIR)/omp_impl.c -o $(BINDIR)/omp_impl.o


# Tests with debug printing 
testd: testd_omp testd_mpi
	# FIXME: How to use bash variables in makefile?
	# for file in "$(TESTFILES)[@]"; do if [ ! -e file ]; then python3 $(TESTDIR)/generate_matrices.py; fi; done

testd_mpi: $(SRCDIR)/test.c $(SRCDIR)/omp_impl.c
	$(CC) $(TESTFLAGS) -c $(SRCDIR)/mpi_impl.c -o $(TESTDIR)/mpi_impl.o
	$(CC) $(TESTFLAGS) -c -D_MPI_TEST $(SRCDIR)/test.c -o $(TESTDIR)/mpi_test.o
	$(MC) $(TESTFLAGS) $(TESTDIR)/mpi_test.o $(TESTDIR)/mpi_impl.o -o $(TESTDIR)/mpi_test.out
	cd $(TESTDIR) && mpiexec mpi_test.out -n 6
	cd ..

testd_omp: $(SRCDIR)/test.c $(SRCDIR)/omp_impl.c
	$(CC) $(TESTFLAGS) $(OMPFLAGS) -c $(SRCDIR)/omp_impl.c -o $(TESTDIR)/omp_impl.o
	$(CC) $(TESTFLAGS) $(OMPFLAGS) -c -D_OMP_TEST $(SRCDIR)/test.c -o $(TESTDIR)/omp_test.o
	$(CC) $(TESTFLAGS) $(OMPFLAGS) $(TESTDIR)/omp_impl.o $(TESTDIR)/omp_test.o -o $(TESTDIR)/omp_test.out
	cd $(TESTDIR) && ./omp_test.out
	cd ..


# Tests without debug printing
test: test_omp test_mpi

test_omp: $(SRCDIR)/omp_impl.c $(SRCDIR)/test.c
	if [ ! -e $(TESTDIR) ]; then mkdir $(TESTDIR); fi
	$(CC) $(CFLAGS) $(OMPFLAGS) $(TESTFLAGS) -DNDEBUG $(SRCDIR)/test.c $(SRCDIR)/omp_impl.c -o $(TESTDIR)/test.out

test_mpi: $(SRCDIR)/mpi_impl.c
	echo "MPI testing not implemented"

clean:
	rm $(BINDIR)/*.out
	rm $(BINDIR)/*.o
	rm $(TESTDIR)/*.out
	rmdir $(BINDIR)
