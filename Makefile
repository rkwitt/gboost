
all:	Makefile.options
	make -C src-graphmatch/ all
	make -C src-gspan/ all

	cp src-graphmatch/graphmatch.mex* bin/
	cp src-gspan/mexgspan.mex* bin/

