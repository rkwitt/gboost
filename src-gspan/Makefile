
###
### Options
###

include ../Makefile.options


###
### The values below should not be modified
###

CPP = g++
INCLUDE=-I${MATLABROOT}/extern/include

# DEBUG and RELEASE settings
#CPPFLAGS = -DDEBUG -g -ggdb -Wall -Wno-deprecated ${INCLUDE}
#CPPFLAGS = -g -ggdb -Wall -Wno-deprecated ${INCLUDE}
CPPFLAGS = -DNDEBUG -O3 -fomit-frame-pointer -funroll-loops -fPIC \
	-Wall -Wno-deprecated ${INCLUDE}

MEX=${MATLABROOT}/bin/mex
MEXFLAGS=-O

OBJ = graph.o misc.o dfs.o ismin.o gspan.o  main.o


### Targets

all:	gspan

gspan:	${OBJ} mexgspan.cpp
	${MEX} ${MEXFLAGS} mexgspan.cpp ${OBJ}

clean:
	rm -f *.o

