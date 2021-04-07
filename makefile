CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic -pg
OMPFLAGS=-fopenmp
LDFLAGS=-lm


SRC=src
OBJ=obj
BIN=bin
TEST=test

SRCS=$(wildcard $(SRC)/*.c)
SRCTEST=$(wildcard $(SRC)/*test.c)
SRCOMP=$(wildcard $(SRC)/*omp.c)
SRCMPI=$(wildcard $(SRC)/*mpi.c)
SRCIMPL=$(wildcard $(SRC)/impl*.c)
OBJOMP=$(patsubst $(SRC)/%omp.c, $(OBJ)/%omp.o, $(SRCOMP))
OBJMPI=$(patsubst $(SRC)/%mpi.c, $(OBJ)/%mpi.o, $(SRCMPI))
OBJIMPL=$(patsubst $(SRC)/impl%.c, $(OBJ)/impl%.o, $(SRCIMPL))

BINOMP=$(BIN)/omp.out
BINMPI=$(BIN)/mpi.out

all:$(BINOMP) # $(BINMPI)


release: CFLAGS=-Wall -pedantic -O2 -DNDEBUG
release: clean
release: $(BINOMP) $(BINMPI)

yeet: $(OBJIMPL) $(SRC)/mpi_test.c $(SRC)/omp_test.c
	$(CC) $(CFLAGS) $(OMPFLAGS) $(SRC)/omp_test.c $(OBJ)/impl_omp.o $(LDFLAGS) -o $(BIN)/test_omp.out
#	$(MC) $(CFLAGS) $(SRC)/test_mpi.c $(OBJ)/impl_mpi.o $(LDFLAGS) -o $(BIN)/test_mpi.out

$(BIN)/test_%.out: $(SRCTEST) $(OBJ)/impl_%.o

$(BINOMP): $(OBJOMP)
	$(CC) $(CFLAGS) $(OMPFLAGS) $(OBJOMP) -o $(BINOMP) $(LDFLAGS)

$(BINMPI): $(OBJMPI)
	# $(MC) $(CFLAGS) $(OBJMPI) -o $(BINMPI)

$(OBJ)/%omp.o: $(SRC)/%omp.c
	$(CC) -c $(CFLAGS) $(OMPFLAGS) $< -o $@ $(LDFLAGS)

$(OBJ)/%mpi.o: $(SRC)/%mpi.c
	#$(MC) -c $(CLAGS) $(OMPFLAGS) $< -o $@

clean:
	rm -r $(BIN)/* $(OBJ)/*

