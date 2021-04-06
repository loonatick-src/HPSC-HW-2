CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic -pg
OMPFLAGS=-fopenmp

SRC=src
OBJ=obj
BIN=bin
TEST=test

SRCS=$(wildcard $(SRC)/*.c)
SRCOMP=$(wildcard $(SRC)/*omp.c)
SRCMPI=$(wildcard $(SRC)/*mpi.c)
OBJOMP=$(patsubst $(SRC)/%omp.c, $(OBJ)/%omp.o, $(SRCOMP))
OBJMPI=$(patsubst $(SRC)/%mpi.c, $(OBJ)/%mpi.o, $(SRCMPI))

BINOMP=$(BIN)/omp.out
BINMPI=$(BIN)/mpi.out

all:$(BINOMP) # $(BINMPI)

release: CFLAGS=-Wall -pedantic -O2 -DNDEBUG
release: clean
release: $(BINOMP) # $(BINMPI)

$(BINOMP): $(OBJOMP)
	$(CC) $(CFLAGS) $(OMPFLAGS) $(OBJOMP) -o $(BINOMP)

$(BINMPI): $(OBJMPI)
	echo "MPI not implemented"	
	# $(MC) $(CFLAGS) $(OBJMPI) -o $(BINMPI)

$(OBJ)/%omp.o: $(SRC)/%omp.c
	$(CC) -c $(CFLAGS) $(OMPFLAGS) $< -o $@

$(OBJ)/%mpi.o: $(SRC)/%mpi.c
	echo "Not Implemented"
	#$(MC) -c $(CLAGS) $(OMPFLAGS) $< -o $@

clean:
	rm -r $(BIN)/* $(OBJ)/*

