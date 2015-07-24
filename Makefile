#
# MatchMaker
#   

INC_FITS=-I/home/dmonge/include
LIB_FITS=-L/home/dmonge/lib
INC_GSL=
LIB_GSL=

CC      = mpicc -std=c99 
WOPT    ?= -Wall
CFLAGS  := -O3 $(WOPT)
#CFLAGS  := -g -O0 $(WOPT)
CFLAGS += -D_LONGIDS -D_DEBUG
CFLAGS += -I./ $(INC_FITS) $(INC_GSL)
LIBS    := $(LIB_FITS) $(LIB_GSL) -lcfitsio -lgsl -lgslcblas -lm

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

clean :
	rm -f $(EXEC) $(OBJS) $(OBJS2)

cleaner :
	rm -f $(EXEC) $(OBJS) $(OBJS2) *~ src/*~
