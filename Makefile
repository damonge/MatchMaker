#
# COLA_DAM
#   

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
#OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc

CC      = mpicc -std=c99 
WOPT    ?= -Wall
CFLAGS  := -O3 $(WOPT) $(OPENMP)
CFLAGS += #-D_DEBUG
LIBS    := -lm

# Define paths of FFTW3 & GSL libraries if necessary.

EXEC = MatchMaker
all: $(EXEC)

OBJS := src/main.o
OBJS += src/read_param.o src/snap_io.o src/msg.o src/fof.o

MatchMaker: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

src/main.o: src/main.c src/common.h
src/read_param.o: src/read_param.c src/common.h
src/snap_io.o: src/snap_io.c src/common.h
src/msg.o: src/msg.c src/common.h
src/fof.o: src/fof.c src/common.h

.PHONY: clean run dependence
clean :
	rm -f $(EXEC) $(OBJS) $(OBJS2)

cleaner :
	rm -f $(EXEC) $(OBJS) $(OBJS2) *~ src/*~

run:
	mpirun -n 2 ./COLA_DAM param.ini

dependence:
	gcc -MM -MG *.c
