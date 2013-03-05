
### Options

CPP=g++
CPPFLAGS=-g -ggdb -Wall

INCLUDE=-I./vflib-2.0/include
LIBS=-L./vflib-2.0/lib -lvf

include ../Makefile.options

MEX=${MATLABROOT}/bin/mex
MEXFLAGS=-O
#MEXFLAGS=-DDEBUG


### Targets

all:	graphmatch.mex

graphmatch.mex:	graphmatch.cpp Makefile
	make -C vflib-2.0 clean all
	$(MEX) $(MEXFLAGS) graphmatch.cpp \
		$(INCLUDE) $(LIBS)

clean:
	make -C vflib-2.0 clean

