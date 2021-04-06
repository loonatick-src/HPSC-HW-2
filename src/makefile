CC=gcc
MC=mpicc
CFLAGS=-Wall -Winline -pedantic -pg
OMPFLAGS=-fopenmp
SRC=src
OBJ=obj
BIN=bin
TEST=test
BINARIES=$(patsubst $($(SRC)/impl%.c), $(BIN)/%.out, )
SRCS=$(wildcard $(SRC)/*.c)
SRCOMP=$(wildcard $(SRC)/*omp.c)
SRCMPI=$(wildcard $(SRC)/*mpi.c)
OBJS=$(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCS))
OBJOMP=$(wildcard $(OBJ)/*omp.o)
OBJMPI=$(wildcard $(OBJ)/*.mpi.o)

all:$(BIN)

$(BIN): 

release: CFLAGS=-Wall -Winline -pedanic -O2 -DNDEBUG
release: clean
release: $(BIN)



clean:
	rm $(OBJS)
	rm $(BIN)/*
