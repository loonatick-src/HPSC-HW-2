CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic  
OMPFLAGS=-fopenmp
DBGFLAGS=-g

omp: main_omp.c omp_impl.o
	$(CC) $(CFLAGS) $(OMPFLAGS) main_omp.c omp_impl.o -o omp.out

omp_impl.o: omp_impl.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c omp_impl.c


